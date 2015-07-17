#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>
#include <math.h>

#include <utils/array2d/carray2d.h>
// #include "utils/vector_utils/vector_utils.h"

// #define STB_IMAGE_IMPLEMENTATION
#include <poro/external/stb_image.h>
// #undef STB_IMAGE_IMPLEMENTATION

// #define STB_IMAGE_WRITE_IMPLEMENTATION
#include <poro/external/stb_image_write.h>
#include "external/jpge.h"

// #undef STB_IMAGE_WRITE_IMPLEMENTATION

//-----------------------------------------------------------------------------
// Image handling stuff

void LoadImageTo( const std::string& filename, ceng::CArray2D< unsigned int >& out_texture_data )
{
	int width,height,bpp;
	unsigned char *data = stbi_load(filename.c_str(), &width, &height, &bpp, 4);
	out_texture_data.Resize( width, height );

	memcpy( out_texture_data.GetData().data, (unsigned int*)data, width * height * sizeof( unsigned int ) );
	stbi_image_free( data );
}

// saves a jpg
void SaveImageTo( const std::string& filename, const ceng::CArray2D< unsigned int >& image_out )
{
	int quality_factor = 70;
	int subsampling = jpge::H2V2;
	bool optimize_huffman_tables = true;

	using namespace jpge;

	const int req_comps = 3; // request RGB image
	int width, height;

	width = image_out.GetWidth();
	height = image_out.GetHeight();

	jpge::params params;
	params.m_quality = quality_factor;
	params.m_subsampling = (subsampling < 0) ? ((req_comps == 1) ? jpge::Y_ONLY : jpge::H2V2) : static_cast<jpge::subsampling_t>(subsampling);
	params.m_two_pass_flag = optimize_huffman_tables;

	unsigned char* image_data = new unsigned char[ width * height * req_comps ];
	unsigned char* imd = image_data;
	unsigned char* source = (unsigned char*)( image_out.GetData().data );


	for( int y = 0; y < height; ++y )
	{
		for( int x = 0; x < width; ++x )
		{
			for( int i = 0; i < 3; ++i )
			{	
				(*imd) = *source;
				++imd;
				++source;
			}
			
			++source;
		}
	}


	if (!jpge::compress_image_to_jpeg_file(filename.c_str(), width, height, req_comps, image_data, params))
	{
		std::cout << "JPGE::Compress_image_to_jpeg_file() failed with filename: " << filename << std::endl;
	}

	delete [] image_data;
}

void SavePngTo( const std::string& filename, ceng::CArray2D< unsigned int >& image_out )
{
	
	const int bpp = 4;

	int result = stbi_write_png( filename.c_str(), image_out.GetWidth(), image_out.GetHeight(), 4, image_out.GetData().data, image_out.GetWidth() * 4 );

}



struct IMPL_Surface {
	unsigned char* pixels;
	int w;
	int h;
	int pitch;
};

typedef struct tColorRGBA {
	unsigned char r;
	unsigned char g;
	unsigned char b;
	unsigned char a;
} tColorRGBA;


int ResizeImage( ceng::CArray2D< unsigned int >& surface, int new_width, int new_height )
{
	ceng::CArray2D< unsigned int > new_surface( new_width, new_height );
	int flipx = 0;
	int flipy = 0;
	int smooth = 1;

	IMPL_Surface surface_src;
	surface_src.w = surface.GetWidth();
	surface_src.h = surface.GetHeight();
	surface_src.pitch = surface_src.w * 4;
	surface_src.pixels = (unsigned char*)( surface.GetData().data );

	IMPL_Surface surface_dst;
	surface_dst.w = new_surface.GetWidth();
	surface_dst.h = new_surface.GetHeight();
	surface_dst.pitch = surface_dst.w * 4;
	surface_dst.pixels = (unsigned char*)( new_surface.GetData().data );

	IMPL_Surface* src = &surface_src;
	IMPL_Surface* dst = &surface_dst;
	
	// int _zoomSurfaceRGB(IMPL_Surface* src, IMPL_Surface* dst, int flipx, int flipy, int smooth)

	const int comp = 4;
	typedef unsigned int Uint32;
	typedef unsigned char Uint8;

	int x, y, sx, sy, ssx, ssy, *sax, *say, *csax, *csay, *salast, csx, csy, ex, ey, cx, cy, sstep, sstepx, sstepy;
	tColorRGBA *c00, *c01, *c10, *c11;
	tColorRGBA *sp, *csp, *dp;
	int spixelgap, spixelw, spixelh, dgap, t1, t2;

	/*
	* Allocate memory for row/column increments 
	*/
	if ((sax = (int *) malloc((dst->w + 1) * sizeof(Uint32))) == NULL) {
		return (-1);
	}
	if ((say = (int *) malloc((dst->h + 1) * sizeof(Uint32))) == NULL) {
		free(sax);
		return (-1);
	}

	/*
	* Precalculate row increments 
	*/
	spixelw = (src->w - 1);
	spixelh = (src->h - 1);
	if (smooth) {
		sx = (int) (65536.0 * (float) spixelw / (float) (dst->w - 1));
		sy = (int) (65536.0 * (float) spixelh / (float) (dst->h - 1));
	} else {
		sx = (int) (65536.0 * (float) (src->w) / (float) (dst->w));
		sy = (int) (65536.0 * (float) (src->h) / (float) (dst->h));
	}

	/* Maximum scaled source size */
	ssx = (src->w << 16) - 1;
	ssy = (src->h << 16) - 1;

	/* Precalculate horizontal row increments */
	csx = 0;
	csax = sax;
	for (x = 0; x <= dst->w; x++) {
		*csax = csx;
		csax++;
		csx += sx;

		/* Guard from overflows */
		if (csx > ssx) { 
			csx = ssx; 
		}
	}

	/* Precalculate vertical row increments */
	csy = 0;
	csay = say;
	for (y = 0; y <= dst->h; y++) {
		*csay = csy;
		csay++;
		csy += sy;

		/* Guard from overflows */
		if (csy > ssy) {
			csy = ssy;
		}
	}

	sp = (tColorRGBA *) src->pixels;
	dp = (tColorRGBA *) dst->pixels;
	dgap = dst->pitch - dst->w * comp;
	spixelgap = src->pitch/comp;

	if (flipx) sp += spixelw;
	if (flipy) sp += (spixelgap * spixelh);

	/*
	* Switch between interpolating and non-interpolating code 
	*/
	if (smooth) {

		/*
		* Interpolating Zoom 
		*/
		csay = say;
		for (y = 0; y < dst->h; y++) {
			csp = sp;
			csax = sax;
			for (x = 0; x < dst->w; x++) {
				/*
				* Setup color source pointers 
				*/
				ex = (*csax & 0xffff);
				ey = (*csay & 0xffff);
				cx = (*csax >> 16);
				cy = (*csay >> 16);
				sstepx = cx < spixelw;
				sstepy = cy < spixelh;
				c00 = sp;
				c01 = sp;
				c10 = sp;
				if (sstepy) {
					if (flipy) {
						c10 -= spixelgap;
					} else {
						c10 += spixelgap;
					}
				}
				c11 = c10;
				if (sstepx) {
					if (flipx) {
						c01--;
						c11--;
					} else {
						c01++;
						c11++;
					}
				}

				/*
				* Draw and interpolate colors 
				*/
				t1 = ((((c01->r - c00->r) * ex) >> 16) + c00->r) & 0xff;
				t2 = ((((c11->r - c10->r) * ex) >> 16) + c10->r) & 0xff;
				dp->r = (((t2 - t1) * ey) >> 16) + t1;
				t1 = ((((c01->g - c00->g) * ex) >> 16) + c00->g) & 0xff;
				t2 = ((((c11->g - c10->g) * ex) >> 16) + c10->g) & 0xff;
				dp->g = (((t2 - t1) * ey) >> 16) + t1;
				t1 = ((((c01->b - c00->b) * ex) >> 16) + c00->b) & 0xff;
				t2 = ((((c11->b - c10->b) * ex) >> 16) + c10->b) & 0xff;
				dp->b = (((t2 - t1) * ey) >> 16) + t1;
				/*t1 = ((((c01->a - c00->a) * ex) >> 16) + c00->a) & 0xff;
				t2 = ((((c11->a - c10->a) * ex) >> 16) + c10->a) & 0xff;
				dp->a = (((t2 - t1) * ey) >> 16) + t1;				*/
				/*
				* Advance source pointer x
				*/
				salast = csax;
				csax++;				
				sstep = (*csax >> 16) - (*salast >> 16);
				if (flipx) {
					sp -= sstep;
				} else {
					sp += sstep;
				}

				/*
				* Advance destination pointer x
				*/
				dp++;
			}
			/*
			* Advance source pointer y
			*/
			salast = csay;
			csay++;
			sstep = (*csay >> 16) - (*salast >> 16);
			sstep *= spixelgap;
			if (flipy) { 
				sp = csp - sstep;
			} else {
				sp = csp + sstep;
			}

			/*
			* Advance destination pointer y
			*/
			dp = (tColorRGBA *) ((Uint8 *) dp + dgap);
		}
	} else {
		/*
		* Non-Interpolating Zoom 
		*/		
		csay = say;
		for (y = 0; y < dst->h; y++) {
			csp = sp;
			csax = sax;
			for (x = 0; x < dst->w; x++) {
				/*
				* Draw 
				*/
				*dp = *sp;

				/*
				* Advance source pointer x
				*/
				salast = csax;
				csax++;				
				sstep = (*csax >> 16) - (*salast >> 16);
				if (flipx) sstep = -sstep;
				sp += sstep;

				/*
				* Advance destination pointer x
				*/
				dp++;
			}
			/*
			* Advance source pointer y
			*/
			salast = csay;
			csay++;
			sstep = (*csay >> 16) - (*salast >> 16);
			sstep *= spixelgap;
			if (flipy) sstep = -sstep;			
			sp = csp + sstep;

			/*
			* Advance destination pointer y
			*/
			dp = (tColorRGBA *) ((Uint8 *) dp + dgap);
		}
	}

	/*
	* Remove temp arrays 
	*/
	free(sax);
	free(say);

	surface = new_surface;
	return (0);
}

//-----------------------------------------------------------------------------


struct Vec2
{
	Vec2() { }
	Vec2( int x, int y ) : x(x), y(y) { }

	bool operator==( const Vec2& other ) const { return ( x == other.x && y == other.y ); }
	bool operator!=( const Vec2& other ) const { return !( *this == other ); }

	void Set( int _x, int _y ) { x = _x; y = _y; }
	int x, y;
};


struct Line
{
	Vec2 a;
	Vec2 b;
};

struct LinePair
{
	Line a;
	Line b;
};

int mod(int x, int m) {
    int r = x%m;
    return r<0 ? r+m : r;
}

template <class T> void swap( T& a, T& b ) {
  T c(a); a=b; b=c;
}

float Length( const Vec2& vec2 ) {
	return sqrtf( (float)( vec2.x * vec2.x + vec2.y * vec2.y ) );
}

float Distance( const Vec2& a, const Vec2& b ) {
	return Length( Vec2( a.x - b.x, a.y - b.y ) );
}

//--------------------------


void WritePng( const ceng::CArray2D< bool >& image, const std::string& filename, std::vector< Vec2 >& positions )
{
	const int bpp = 4;
	unsigned char* data = new unsigned char[ image.GetWidth() * image.GetHeight() * 4 ];

	for( int y = 0; y < image.GetHeight(); ++y )
	{
		for( int x = 0; x < image.GetWidth(); ++x )
		{
			bool is_marker = false;
			if( std::find( positions.begin(), positions.end(), Vec2( x, y ) ) != positions.end() )
				is_marker = true;

			int i = ( y * image.GetWidth() + x ) * bpp;
			unsigned char value = image.At( x, y ) ? 255 : 0;
			data[i] = value;
			data[i+1] = value;
			data[i+2] = value;
			data[i+3] = value;

			if( is_marker ) 
			{
				data[i] = 255;
				data[i+1] = 128;
				data[i+2] = 0;
				data[i+3] = 255;
			}
		}
	}

	int result = stbi_write_png( filename.c_str(), image.GetWidth(), image.GetHeight(), 4, data, image.GetWidth() * 4 );

	delete [] data;
}

void ParseBlackAndWhite( unsigned char* image, int w, int h, int bpp, unsigned char threshold, ceng::CArray2D< bool >& output )
{
	output.Resize( w, h );
	output.SetEverythingTo( false );
	
	for( int y = 0; y < h; ++y )
	{
		for( int x = 0; x < w; ++x )
		{
			int i = ( y * w + x ) * bpp;
			if( ( image[i] > threshold ) || ( image[i+1] > threshold ) || ( image[i+2] > threshold ) )
				output.At( x, y ) = true;
		}
	}
}

// 95% percent compare
bool CmprN( int n1, int n2 )
{
	if( n1 == 0 || n2 == 0 ) return false;
	if( n1 == n2 ) return true;
	if( n1 > n2 ) swap( n1, n2 );
	const int diff = n2 - n1;
	const int max_allowed = (int) ( ((float)(n1) * 0.90f) );
	if( diff <= max_allowed ) return true;
	return false;
}

int FindNextLine( const std::vector< Line >& lines, const std::vector< bool >& used, const Vec2& start, const Vec2& end, int i )
{
	float delta_x = (float)( end.x - start.x );
	float delta_y = (float)( end.y - start.y );
	float length = Length( Vec2( end.x - start.x, end.y - start.y ) );
	if( length != 0 )
	{
		delta_x /= length;
		delta_y /= length;
	}

	const int limit = 3;
	const float delta_limit = 1.75f;

	float best_delta = 1.75f;
	float closest = limit;
	int best_i = -1;

	for( ; i < (int)lines.size(); ++i )
	{
		if( used[i] ) continue;
		Vec2 t = lines[i].a;
		if( abs( end.x - t.x ) > limit ||
			abs( end.y - t.y ) > limit ) continue;

		// check the delta
		// float t_delta_x = (float)( t.x - end.x );
		// float t_delta_y = (float)( t.y - end.y );
		float t_length = Length( Vec2( t.x - end.x, t.y - end.y ) );
		// t_delta_x /= t_length;
		// t_delta_y /= t_length;
		
		/*float diff = abs( t_delta_x - delta_x );
		if( abs( t_delta_y - delta_y ) > diff ) 
			diff = abs( t_delta_y - delta_y );*/

		float distance = t_length;
		if( distance < closest )
		{
			closest = distance;
			// best_delta = diff;
			best_i = i;
		}
	}

	if( length == 0 || closest < limit ) 
		return best_i;

	return -1;
}

std::vector< LinePair > mLines;

void FindLines( const std::vector< Line >& lines )
{
	std::vector< bool > used( lines.size() );
	std::vector< LinePair >& results = mLines;
	for( std::size_t i = 0; i < lines.size(); ++i )
	{
		used[i] = true;
		Line start = lines[i];
		Line end = lines[i];
		int next = FindNextLine( lines, used, start.a, end.a, i + 1 );
		while( next != -1 )
		{
			end = lines[ next ];
			used[ next ] = true;
			next = FindNextLine( lines, used, start.a, end.a, i + 1 );
		}
		
		if( start.a != end.a && start.b != end.b )
		{
			LinePair l;
			l.a.a = start.a;
			l.a.b = end.a;
			l.b.a = start.b;
			l.b.b = end.b;
			results.push_back( l );
		}
	}

	std::cout << "Found lines: " << results.size() << std::endl;
}


void LookForMarkers( const ceng::CArray2D< bool >& image, std::vector< Vec2 >& out_positions )
{

	std::vector< Line > horizontal;
	std::vector< Line > vertical;
	// 0,  1,  0,  1,  0 
	// 1 : 1 : 3 : 1 : 1
	// ring buffer
	std::vector< Vec2 > start_positions( 6 );
	std::vector< Vec2 > end_positions( 6 );
	std::vector< int > lengths( 6 );
	int li = 0;
	bool looking_for = false;
	start_positions[0].Set( 0, 0 );

	for( int y = 0; y < image.GetHeight(); ++y )
	{
		memset( &lengths[0], 0, lengths.size() * sizeof( int ) );
		li = 0;
		looking_for = false;
		start_positions[0].Set( 0, y );

		for( int x = 0; x < image.GetWidth(); ++x )
		{
			bool v = image.At( x, y );
			if( v == looking_for )
			{
				lengths[li]++;
			}
			else
			{
				end_positions[li].Set( x - 1, y );
				li++;
				if( li >= (int)lengths.size() ) li = 0;
				lengths[li] = 0;
				looking_for = !looking_for;

				start_positions[li].Set( x, y );

				// check if we have a markers?
				int l1 = lengths[ mod( li - 1, (int)lengths.size() ) ];
				int l2 = lengths[ mod( li - 2, (int)lengths.size() ) ];
				int l3 = lengths[ mod( li - 3, (int)lengths.size() ) ];
				int l4 = lengths[ mod( li - 4, (int)lengths.size() ) ];
				int l5 = lengths[ mod( li - 5, (int)lengths.size() ) ];

				if( CmprN( l1, l2 ) && CmprN( l1 * 3, l3 ) && CmprN( l2, l4 ) && CmprN( l4, l5 ) )
				{
					// found a marker
					// std::cout << "found a marker at "  << x << ", " << y << std::endl;
					
					out_positions.push_back( start_positions[ mod( li - 5, (int)start_positions.size() ) ] );
					out_positions.push_back( Vec2( x - 1, y ) );

					Line t;
					t.a = start_positions[ mod( li - 5, (int)start_positions.size() ) ];
					t.b = Vec2( x - 1, y );
					horizontal.push_back( t );
				}
			}
		}
	}

	// --- copy of above --
	for( int x = 0; x < image.GetWidth(); ++x )
	{
		memset( &lengths[0], 0, lengths.size() * sizeof( int ) );
		li = 0;
		looking_for = false;
		start_positions[0].Set( x, 0 );

		for( int y = 0; y < image.GetHeight(); ++y )
		{
			bool v = image.At( x, y );
			if( v == looking_for )
			{
				lengths[li]++;
			}
			else
			{
				end_positions[li].Set( x, y - 1 );
				li++;
				if( li >= (int)lengths.size() ) li = 0;
				lengths[li] = 0;
				looking_for = !looking_for;

				start_positions[li].Set( x, y );

				// check if we have a markers?
				int l1 = lengths[ mod( li - 1, (int)lengths.size() ) ];
				int l2 = lengths[ mod( li - 2, (int)lengths.size() ) ];
				int l3 = lengths[ mod( li - 3, (int)lengths.size() ) ];
				int l4 = lengths[ mod( li - 4, (int)lengths.size() ) ];
				int l5 = lengths[ mod( li - 5, (int)lengths.size() ) ];

				if( CmprN( l1, l2 ) && CmprN( l1 * 3, l3 ) && CmprN( l2, l4 ) && CmprN( l4, l5 ) )
				{
					// found a marker
					// std::cout << "found a marker at "  << x << ", " << y << std::endl;
					
					out_positions.push_back( start_positions[ mod( li - 5, (int)start_positions.size() ) ] );
					out_positions.push_back( Vec2( x, y - 1 ) );

					Line t;
					t.a = start_positions[ mod( li - 5, (int)start_positions.size() ) ];
					t.b = Vec2( x, y - 1 );
					vertical.push_back( t );
				}
			}
		}
	}

	// Look for points that are next to each other
	std::vector< Line >& lines = horizontal;
	std::cout << "points: " << horizontal.size() << std::endl;
	FindLines( horizontal );
	FindLines( vertical );

}

ceng::CArray2D< bool > image;
unsigned int* image_data = NULL;
ceng::CArray2D< unsigned int > mImage;


int ParseQR_LoadImage( const std::string& filename )
{
	LoadImageTo( filename, mImage );
	int work_width = 1500;

	if( mImage.GetWidth() == 0 ) 
	{
		std::cout << "Error couldn't load file: " << filename << std::endl;
		return -1;
	}

	if( mImage.GetWidth() != work_width )
	{
		int height = (int)( (float)( (float)work_width / (float)mImage.GetWidth() ) * (float)mImage.GetHeight() );
		ResizeImage( mImage, work_width, height );
	}

	for( int y = 0; y < mImage.GetHeight(); ++y )
	{
		for( int x = 0; x < mImage.GetWidth(); ++x )
		{
			mImage.At( x, y ) = ( mImage.At( x, y ) & 0x00FFFFFF ) | 0xFF000000;
		}
	}

	image_data = mImage.GetData().data;

	return 0;
}

int ParseQR( const std::string& filename, int iteration )
{
	// Step 1: Parse a black and white image 
	// Step 2: Look for markers

	if( image_data == NULL )
	{
		ParseQR_LoadImage( filename );
	}

	/*
	int x,y,bpp;
	unsigned char *data = stbi_load(filename.c_str(), &x, &y, &bpp, 4);
	image_data = (unsigned int*)data;

	if( data == NULL ) {
		std::cout << "Error! - Couldn't read file: " << filename << std::endl;
		return -1;
	}
	*/

	// check for possible places
	unsigned char thresholds[] = { 128, 
		128 + 64, 128 - 64,
		128 + 32, 128 - 32, 128 + 64 + 32, 128 - 64 - 32,
		128 + 16, 128 - 16, 128 + 64 + 16, 128 - 64 - 16, 128 + 64 + 32 + 16, 128 - 64 - 32 - 16 };


	unsigned char threshold = thresholds[ iteration ];
	// threshold = 96;
	std::cout << "testing treshold: " << (int)threshold << std::endl;
	ParseBlackAndWhite( (unsigned char*)image_data, mImage.GetWidth(), mImage.GetHeight(), 4, threshold, image );

	std::vector< Vec2 > positions;
	LookForMarkers( image, positions );

	std::stringstream ss;
	ss << "out_black_n_white_" << (int)threshold << ".png";
	WritePng( image, ss.str(), positions );

	
	return 1;
}

#if 0

int main( int argc, char** args )
{
	ParseQR( "test/test_006.png" );
	return 0;
}

#endif