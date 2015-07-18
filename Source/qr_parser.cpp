#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>
#include <math.h>
#include <float.h>

#define COMMAND_LINE
// #define QR_DEBUG
#include "external/jpge.h"

#ifdef COMMAND_LINE
	#include "qr_utils.h"
	
	#define STB_IMAGE_IMPLEMENTATION
	#include "external/stb_image.h"
	#undef STB_IMAGE_IMPLEMENTATION

	#include "external/jpge.cpp"

	// #define QR_DEBUG
	#ifdef QR_DEBUG
	#define STB_IMAGE_WRITE_IMPLEMENTATION
	#include "external/stb_image_write.h"
	#endif

#endif


#ifndef COMMAND_LINE
	// TestAABBIntersection
	
	#include <utils/math/math_utils.h>
	#include <utils/math/point_inside.h>
	#include <utils/color/ccolor.h>
	#include <utils/random/random.h>
	#include <utils/array2d/carray2d.h>
	#include <utils/staticarray/cstaticarray.h>

	// #define STB_IMAGE_IMPLEMENTATION
	#include <poro/external/stb_image.h>
	// #undef STB_IMAGE_IMPLEMENTATION

	// #define QR_DEBUG
	#ifdef QR_DEBUG
	#include <poro/external/stb_image_write.h>
	#endif

#endif 


// #undef STB_IMAGE_WRITE_IMPLEMENTATION

//-----------------------------------------------------------------------------
// Image handling stuff

void LoadImageTo( const std::string& filename, ceng::CArray2D< unsigned int >& out_texture_data )
{
	int width,height,bpp;
	unsigned char *data = stbi_load(filename.c_str(), &width, &height, &bpp, 4);
	
	if( data == NULL || width <= 0 || height <= 0 )
	{
		std::cout << "Error couldn't load file: " << filename << std::endl;
	}
	else
	{
		out_texture_data.Resize( width, height );

		memcpy( out_texture_data.GetData().data, (unsigned int*)data, width * height * sizeof( unsigned int ) );
		stbi_image_free( data );
	}
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

#ifdef QR_DEBUG
void SavePngTo( const std::string& filename, ceng::CArray2D< unsigned int >& image_out )
{
	const int bpp = 4;
	int result = stbi_write_png( filename.c_str(), image_out.GetWidth(), image_out.GetHeight(), 4, image_out.GetData().data, image_out.GetWidth() * 4 );
}
#endif



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
	if( surface.GetWidth() <= 0 || surface.GetHeight() <= 0 ||
		new_width <= 0 || new_height <= 0 )
	{
		return -1;
	}

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
	Vec2() : x(0), y(0) { }
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


enum CORNERS
{
	C_TOP_LEFT = 0,
	C_TOP_RIGHT = 1,
	C_BOTTOM_RIGHT = 2,
	C_BOTTOM_LEFT = 3,
};

struct Marker
{
	ceng::CStaticArray< types::vector2, 4 > corners;
	types::vector2 center;

	types::vector2 aabb_min;
	types::vector2 aabb_max;

	void UpdateAABB()
	{
		aabb_min.Set( FLT_MAX, FLT_MAX );
		aabb_max.Set( FLT_MIN, FLT_MIN );
		for( int i = 0; i < corners.length; ++i )
		{
			aabb_min.x = ceng::math::Min( aabb_min.x, corners[i].x );
			aabb_min.y = ceng::math::Min( aabb_min.y, corners[i].y );
			aabb_max.x = ceng::math::Max( aabb_max.x, corners[i].x );
			aabb_max.y = ceng::math::Max( aabb_max.y, corners[i].y );
		}
	}
};

//-----------------------------------------------------------------------------

std::vector< LinePair > mLines;
ceng::CArray2D< bool > image;
unsigned int* image_data = NULL;
ceng::CArray2D< unsigned int > mImage;
std::vector< Marker > mMarkers;
std::vector< Marker > mDebugMarkers;

//--------------------------

#ifdef QR_DEBUG
void WritePng( const ceng::CArray2D< bool >& image, const std::string& filename, std::vector< Vec2 >& positions )
{
	const int bpp = 4;
	unsigned char* data = new unsigned char[ image.GetWidth() * image.GetHeight() * 4 ];

	for( int y = 0; y < image.GetHeight(); ++y )
	{
		for( int x = 0; x < image.GetWidth(); ++x )
		{
			bool is_marker = false;
			/*if( std::find( positions.begin(), positions.end(), Vec2( x, y ) ) != positions.end() )
				is_marker = true;*/

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
	
	unsigned int* data_32 = (unsigned int*)data;
	for( std::size_t i = 0; i < positions.size(); ++i )
	{
		data_32[ positions[i].y * image.GetWidth() + positions[i].x ] = 0xFF9900FF;
	}

	int result = stbi_write_png( filename.c_str(), image.GetWidth(), image.GetHeight(), 4, data, image.GetWidth() * 4 );

	delete [] data;
}
#endif // QR_DEBUG

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
	const int max_allowed = (int) ( ((float)(n1) * 0.85f) );
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


void FindLines( const std::vector< Line >& lines )
{
	std::vector< bool > used( lines.size() );
	for( std::size_t i = 0; i < used.size(); ++i ) {
		used[i] = false;
	}

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

				int total_l = l1+l2+l3+l4+l5;

				if( CmprN( l1, l5 ) && CmprN( l2, l4 ) && l3 > l1 && l3 < (l1+l2+l4+l5) && looking_for )
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

				int total_l = l1+l2+l3+l4+l5;

				if( CmprN( l1, l5 ) && CmprN( l2, l4 ) && l3 > l1 && l3 < (l1+l2+l4+l5) && looking_for )
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



int ParseQR_LoadImage( const std::string& filename )
{
	LoadImageTo( filename, mImage );
	// int work_width = 1555;
	int work_width = 1481;

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

#ifdef QR_DEBUG
	SavePngTo( "debug_resized.png", mImage ); 
#endif

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
		128 + 16, 128 - 16, 128 + 64 + 16, 128 - 64 - 16, 128 + 64 + 32 + 16, 128 - 64 - 32 - 16,
		128 + 64 + 32 + 16 + 8, 128 - 64 - 32 - 16 - 8, 128 + 64 + 32 + 16 + 8 + 4, 128 - 64 - 32 - 16 - 8 - 4 };


	unsigned char threshold = thresholds[ iteration ];
	// threshold = 96;
	std::cout << "testing treshold: " << (int)threshold << std::endl;
	ParseBlackAndWhite( (unsigned char*)image_data, mImage.GetWidth(), mImage.GetHeight(), 4, threshold, image );

	std::vector< Vec2 > positions;
	LookForMarkers( image, positions );

#ifdef QR_DEBUG

	std::stringstream ss;
	ss << "out_black_n_white_" << (int)threshold << ".png";
	WritePng( image, ss.str(), positions );
#endif
	
	return 1;
}

void RemoveShortest( std::vector< LinePair >& line_pairs )
{
	if( line_pairs.empty() ) return;

	float shortest = FLT_MAX;
	int shortest_i = -1;
	for( int i = 0; i < (int)line_pairs.size(); ++i )
	{
		float length = Distance( line_pairs[i].a.a, line_pairs[i].a.b );
		if( length < shortest )
		{
			shortest_i = i;
			shortest = length;
		}
	}

	if( shortest_i != -1 )
	{
		line_pairs[shortest_i] = line_pairs[ line_pairs.size() - 1 ];
		line_pairs.pop_back();
	}
}


types::vector2 FindClosestTo( const ceng::CStaticArray< types::vector2, 4 >& points, const types::vector2& to_what )
{
	float min_dist = FLT_MAX;
	int min_i = -1;

	for( int i = 0; i < points.length; ++i )
	{
		float distance = ( points[i] - to_what ).LengthSquared();
		if( distance < min_dist ) 
		{
			min_i = i;
			min_dist = distance;
		}
	}

	if( min_i != -1 )
		return points[min_i];

	return types::vector2(0,0);
}

void ExpandCorner( Marker& marker, int corner, float expansion_rate )
{

	int prev_i = mod( corner - 1, marker.corners.length );
	int next_i = mod( corner + 1, marker.corners.length );

	types::vector2 p = marker.corners[ corner ];
	types::vector2 delta_0 = p - marker.corners[ prev_i ];
	float length_0 = delta_0.Length();
	delta_0 = delta_0.Normalize();

	types::vector2 delta_1 = p - marker.corners[ next_i ];
	float length_1 = delta_1.Length();
	delta_1 = delta_1.Normalize();

	delta_0 = delta_0 * ( length_0 * expansion_rate );
	p += delta_0;

	delta_1 = delta_1 * ( length_1 * expansion_rate );
	p += delta_1;

	marker.corners[ corner ] = p;
	marker.UpdateAABB();
}

void SortOutMarkerCorners( Marker& marker )
{
	ceng::CStaticArray< types::vector2, 4 > new_corners;
	marker.UpdateAABB();
	types::vector2 aabb_min = marker.aabb_min;
	types::vector2 aabb_max = marker.aabb_max;

	new_corners[ C_TOP_LEFT ] = FindClosestTo( marker.corners, types::vector2( aabb_min.x, aabb_min.y ) );
	new_corners[ C_TOP_RIGHT ] = FindClosestTo( marker.corners, types::vector2( aabb_max.x, aabb_min.y ) );
	new_corners[ C_BOTTOM_RIGHT ] = FindClosestTo( marker.corners, types::vector2( aabb_max.x, aabb_max.y ) );
	new_corners[ C_BOTTOM_LEFT ] = FindClosestTo( marker.corners, types::vector2( aabb_min.x, aabb_max.y ) );

	marker.corners = new_corners;
}


void ExtendLine( types::vector2& a, types::vector2& b, float length )
{
	types::vector2 delta = b - a;
	delta = delta.Normalize();
	b += delta * length;
	a -= delta * length;
}


bool AddAsMarker( const LinePair& l0, const LinePair& l1 )
{
	using namespace types;
	using namespace ceng::math;

	Marker result;

	vector2 top_start( l0.a.a );
	vector2 top_end( l0.a.b );
	ExtendLine( top_start, top_end, 350.f );

	vector2 bottom_start( l0.b.a );
	vector2 bottom_end( l0.b.b );
	ExtendLine( bottom_start, bottom_end, 350.f );

	vector2 left_start( l1.a.a );
	vector2 left_end( l1.a.b );
	ExtendLine( left_start, left_end, 350.f );

	vector2 right_start( l1.b.a );
	vector2 right_end( l1.b.b );
	ExtendLine( right_start, right_end, 350.f );

	types::vector2 top_center = Lerp( top_start, top_end, 0.5f );
	types::vector2 bottom_center = Lerp( bottom_start, bottom_end, 0.5f );
	types::vector2 left_center = Lerp( left_start, left_end, 0.5f );
	types::vector2 right_center = Lerp( right_start, right_end, 0.5f );

	// look for the intersection
	if( LineIntersection( top_start, top_end, left_start, left_end, result.corners[C_TOP_LEFT] ) == false )
	{
		// std::cout << "Error couldn't find top left for marker..." << std::endl;
		return false;
	}

	if( LineIntersection( top_start, top_end, right_start, right_end, result.corners[C_TOP_RIGHT] ) == false )
	{
		// std::cout << "Error couldn't find top left for marker..." << std::endl;
		return false;
	}

	if( LineIntersection( bottom_start, bottom_end, right_start, right_end, result.corners[C_BOTTOM_RIGHT] ) == false )
	{
		// std::cout << "Error couldn't find top left for marker..." << std::endl;
		return false;
	}

	if( LineIntersection( bottom_start, bottom_end, left_start, left_end, result.corners[C_BOTTOM_LEFT] ) == false )
	{
		// std::cout << "Error couldn't find top left for marker..." << std::endl;
		return false;
	}
	
	if( LineIntersection( top_center, bottom_center, left_center, right_center, result.center ) == false )
	{
		// std::cout << "Error couldn't find center for marker" << std::endl;
		return false;
	}

	// sort these out...
	SortOutMarkerCorners( result );

	mMarkers.push_back( result );
	return true;
}

const Marker& FindClosestTo( const std::vector< Marker >& markers, const types::vector2& pos )
{
	float min_dist = FLT_MAX;
	int min_i = -1;
	for( std::size_t i = 0; i < markers.size(); ++i )
	{
		float distance = ( pos - markers[i].center ).LengthSquared();
		if( distance < min_dist )
		{
			min_i = (int)i;
			min_dist = distance;
		}
	}

	if( min_i != -1 ) 
		return markers[min_i];

	return markers[0];
}

void FillTheBlanks( ceng::CArray2D< unsigned int >& texture, Marker& marker )
{
	// aabb
	marker.UpdateAABB();
	types::vector2 aabb_min = marker.aabb_min;
	types::vector2 aabb_max = marker.aabb_max;

	aabb_min -= types::vector2( 4, 4 );
	aabb_max += types::vector2( 4, 4 );

	if( aabb_min.x < 0 ) aabb_min.x = 0;
	if( aabb_min.y < 0 ) aabb_min.y = 0;
	if( (int)aabb_max.x >= texture.GetWidth() ) aabb_max.x = (float)texture.GetWidth() - 1;
	if( (int)aabb_max.y >= texture.GetHeight() ) aabb_max.y = (float)texture.GetHeight() - 1;

	const int mask_size = 4;

	unsigned int c = 0xFFFFFFFF;
	for( int y = (int)aabb_min.y; y <= (int)aabb_max.y; ++y )
	{
		for( int x = (int)aabb_min.x; x <= (int)aabb_max.x; ++x )
		{
			if( image.IsValid( x, y ) )
			{
				bool color_me = false;
				if( image.At( x, y ) == false )
					color_me = true;

				for( int iy = -mask_size; iy <= mask_size; ++iy )
				{
					for( int ix = -mask_size; ix <= mask_size; ++ix )
					{
						if( image.IsValid( x + ix, y + iy ) && image.At( x + ix, y + iy ) == false )
						{
							color_me = true;
							goto GOTOLABEL_FINISHED;
						}
					}
				}

GOTOLABEL_FINISHED:;

				if( color_me )
				{
					if( true )
					{
						float ix = (float)( ( x - aabb_min.x ) / ( aabb_max.x - aabb_min.x ) );
						float iy = (float)( ( y - aabb_min.y ) / ( aabb_max.y - aabb_min.y ) );

						if( ceng::Randomf(0.f, 1.f) < 0.5f )
						{
							int pos_x = x + ceng::Random( -5, 5 );
							int pos_y = (int)( ( ceng::Randomf( 0.f, 1.f ) > iy ) ? aabb_min.y - 2 : aabb_max.y + 2 );
							if( pos_x >= 0 && pos_x < image.GetWidth() && pos_y >= 0 && pos_y < image.GetHeight() )
								c = image_data[ pos_y * image.GetWidth() + pos_x ];
						}
						else
						{
							int pos_x = int( ( ceng::Randomf( 0.f, 1.f ) > ix ) ? aabb_min.x - 2 : aabb_max.x + 2 );
							int pos_y = y + ceng::Random( -5, 5 );
							if( pos_x >= 0 && pos_x < image.GetWidth() && pos_y >= 0 && pos_y < image.GetHeight() )
								c = image_data[ pos_y * image.GetWidth() + pos_x ];
						}
					}

					// find the closest edge?
					texture.At( x, y ) = c;
				}
			}
		}
	}
}


void GaussBlur( float t, ceng::CArray2D< unsigned int >& imagedata, int aabb_min_x, int aabb_min_y, int aabb_max_x, int aabb_max_y )
{
	const int GAUSSIAN_RANGE = 2;
	float radius = t;

	int mask_center = (int)( radius * GAUSSIAN_RANGE ); 
	int fmask_size = 2 * mask_center + 1; 
	std::vector< float > fmask( fmask_size );

	// calculate the matrix
	float total = 0;
	for( int i = 0; i <= mask_center; ++i ) 
	{ 
		float val = exp( - i * i / ( 2 * radius * radius ) ); 
		if( i ) 
			total += val; 
		total += val; 
		fmask[ mask_center + i ] = val;
		fmask[ mask_center - i ] = val; 
	}

	for( size_t i = 0; i < fmask.size(); ++i )
		fmask[ i ] /= total;

	ceng::CArray2D< unsigned int > temp_image( imagedata.GetWidth(), imagedata.GetHeight() );

	/*
	for( int y = 0; y < imagedata.GetHeight(); ++y )
	{
		for( int x = 0; x < imagedata.GetWidth(); ++x )
		{
			if(  imagedata.At( x, y ) == 0 )
			{
				types::fcolor color;
				color.Set32( image_data[ y * imagedata.GetWidth() + x ] );
				color.SetA( 0 );
				imagedata.At( x, y ) = color.Get32();
			}
		}
	}*/

	if( aabb_min_x < 0 ) aabb_min_x = 0;
	if( aabb_min_y < 0 ) aabb_min_y = 0;
	if( aabb_max_x >= imagedata.GetWidth() ) aabb_max_x = imagedata.GetWidth() - 1;
	if( aabb_max_y >= imagedata.GetHeight() ) aabb_max_y = imagedata.GetHeight() - 1;

	// horizontal pass
	for( int y = (int)aabb_min_y; y <= (int)aabb_max_y; ++y )
	{
		for( int x = (int)aabb_min_x; x <= (int)aabb_max_x; ++x )
		{
			// temp_image[ x ][ y ] = 0;
			types::fcolor color( 0,0,0,0 );

			/*
			if(  imagedata.At( x, y ) == 0 )
			{
				color.Set32( image_data[ y * imagedata.GetWidth() + x ] );
				color.SetA( 0 );
			}*/

			for( int i = 0; i < (int)fmask.size(); ++i )
			{
				color += ( 
					types::fcolor( imagedata.At( 
						ceng::math::Clamp( x - mask_center + i, 0, temp_image.GetWidth() ), y ) ) * fmask[ i ]  );
			}
			temp_image.At( x, y ) = color.Get32();
		}
	}

	// vertical pass
	for( int y = (int)aabb_min_y; y <= (int)aabb_max_y; ++y )
	{
		for( int x = (int)aabb_min_x; x <= (int)aabb_max_x; ++x )
		{
			types::fcolor color( 0,0,0,0 );
			
			/*
			if( ( temp_image.At( x, y ) & ~(types::fcolor::AMask) ) == 0 )
			{
				color.Set32( image_data[ y * imagedata.GetWidth() + x ] );
				color.SetA( 0 );
			}*/


			for( int i = 0; i < (int)fmask.size(); ++i )
			{
				color += ( 
					types::fcolor(
						temp_image.At( x, ceng::math::Clamp( y - mask_center + i, 0, temp_image.GetHeight() ) ) ) * fmask[ i ] );
			}

			temp_image[ x ][ y ] = color.Get32();
		}
	} 

	// copy the temp image to the real image
	imagedata = temp_image;
}


void Blur( float t, ceng::CArray2D< unsigned int >& imagedata, const types::vector2& aabb_min, const types::vector2& aabb_max )
{
	const int aabb_min_x = (int)aabb_min.x;
	const int aabb_min_y = (int)aabb_min.y;
	const int aabb_max_x = (int)aabb_max.x;
	const int aabb_max_y = (int)aabb_max.y;

	for( int y = aabb_min_y; y < aabb_max_y; ++y )
	{
		for( int x = aabb_min_x; x < aabb_max_x; ++x )
		{

			
			// temp_image[ x ][ y ] = 0;
			types::fcolor color( 0,0,0,0 );
			color += 1.5f * types::fcolor( imagedata.At( x, y ) );
			color += types::fcolor( imagedata.At( x - 1,  y ) );
			color += types::fcolor( imagedata.At( x + 1,  y ) );
			color += types::fcolor( imagedata.At( x,  y - 1 ) );
			color += types::fcolor( imagedata.At( x,  y + 1 ) );

			color *= 1 / 5.5f;

			imagedata.At( x, y ) = color.Get32();
			/*
			if(  imagedata.At( x, y ) == 0 )
			{
				color.Set32( image_data[ y * imagedata.GetWidth() + x ] );
				color.SetA( 0 );
			}*/

		}
	}
}

float sign( const types::vector2& p1, const types::vector2& p2, const types::vector2& p3) {
    return (p1.x - p3.x) * (p2.y - p3.y) - (p2.x - p3.x) * (p1.y - p3.y);
}

float Dot(const types::vector2& a, const types::vector2& b) {
	return a.x * b.x + a.y * b.y;
}


bool PointInTriangle(const types::vector2& pt, const types::vector2& v1, const types::vector2& v2, const types::vector2& v3)
{
    bool b1, b2, b3;

    b1 = sign(pt, v1, v2) < 0.0f;
    b2 = sign(pt, v2, v3) < 0.0f;
    b3 = sign(pt, v3, v1) < 0.0f;

    return ((b1 == b2) && (b2 == b3));
}

types::fcolor SafeAt( const ceng::CArray2D< unsigned int >& input_image, int x, int y )
{
	if( input_image.IsValid( x, y ) == false ) 
		return types::fcolor(1.f,1.f,1.f,1.f);
	return types::fcolor( input_image.At( x, y ) );
}

unsigned int ReadTexture( const ceng::CArray2D< unsigned int >& input_image, 
				  const std::vector< types::vector2 >& input_pos,
				  float u, float v, float w )
{
	using namespace ceng::math;
	types::vector2 p = u * input_pos[0] + v * input_pos[1] + w * input_pos[2];
	
		
	int px = (int)floor(p.x);
	int py = (int)floor(p.y);
	float fx = (float)( p.x - (double)px ); 
	float fy = (float)( p.y - (double)py );
	
	/*fx = fx*fx*(3.0f-2.0f*fx);	
	fy = fy*fy*(3.0f-2.0f*fy);	*/
	
	types::fcolor color = Lerp(
					Lerp( SafeAt( input_image, px, py     ),	SafeAt( input_image, px + 1, py     ), fx ),
					Lerp( SafeAt( input_image, px, py + 1 ),	SafeAt( input_image, px + 1, py + 1 ), fx ), fy );

	return color.Get32();
}


void DrawTriangle( const ceng::CArray2D< unsigned int >& input_image, 
				  const std::vector< types::vector2 >& input_pos,
				  const std::vector< types::vector2 >& output_pos,
				  ceng::CArray2D< unsigned int >& output_image )
{
	using namespace types;
	vector2 aabb_min( FLT_MAX, FLT_MAX );
	vector2 aabb_max( FLT_MIN, FLT_MIN );

	for( std::size_t i = 0; i < input_pos.size(); ++i )
	{
		aabb_min.x = ceng::math::Min( aabb_min.x, output_pos[i].x );
		aabb_min.y = ceng::math::Min( aabb_min.y, output_pos[i].y );
		aabb_max.x = ceng::math::Max( aabb_max.x, output_pos[i].x );
		aabb_max.y = ceng::math::Max( aabb_max.y, output_pos[i].y );
	}

	// barycentric
	const vector2 v0 = output_pos[1] - output_pos[0];
	const vector2 v1 = output_pos[2] - output_pos[0];
	const float d00 = Dot(v0, v0);
	const float d01 = Dot(v0, v1);
	const float d11 = Dot(v1, v1);
	const float denom = d00 *  d11 - d01 * d01;

	for( int y = (int)aabb_min.y; y <= (int)aabb_max.y; ++y )
	{
		for( int x = (int)aabb_min.x; x <= (int)aabb_max.x; ++x )
		{
			const vector2 p( (float)x, (float)y );
			if( PointInTriangle(
				p, 
				output_pos[0], output_pos[1], output_pos[2] ) )
			{
				// we're inside
				const vector2 v2 = p - output_pos[0];
				const float d20 = Dot( v2, v0 );
				const float d21 = Dot( v2, v1 );
				
				float v = (d11 * d20 - d01 * d21) / denom;
				float w = (d00 * d21 - d01 * d20) / denom;
				float u = 1.f - v - w;

				// output_image.At( x, y ) = input_image.At( x % input_image.GetWidth(), y % input_image.GetHeight() );
				output_image.At( x, y ) = ReadTexture( input_image, input_pos, u, v, w );
			}
		}
	}
}

void BlitTo( const ceng::CArray2D< unsigned int >& image_source, ceng::CArray2D< unsigned int >& dest, int pos_x1, int pos_y1, int pos_x2, int pos_y2 )
{
	for( int y = pos_y1; y <= pos_y2; ++y )
	{
		for( int x = pos_x1; x <= pos_x2; ++x )
		{
			types::fcolor color( image_source.At( x, y ) );
			if( color.GetA() <= 0.f )
				continue;

			color.SetA( ceng::math::Clamp( color.GetA() + 0.25f, 0.f, 1.f ) );

			if( color.GetA() >= 1.f )
			{
				dest.At( x, y ) = image_source.At( x, y );
			}
			else
			{
				types::fcolor result = ceng::math::Lerp( 
					types::fcolor( dest.At( x, y ) ), 
					color, color.GetA() );
				result.SetA( 1.f );
				dest.At( x, y ) = result.Get32();
			}
		}
	}
}


void BlitMultiply( const ceng::CArray2D< unsigned int >& image_source, ceng::CArray2D< unsigned int >& dest, int pos_x1, int pos_y1, int pos_x2, int pos_y2 )
{
	for( int y = pos_y1; y <= pos_y2; ++y )
	{
		for( int x = pos_x1; x <= pos_x2; ++x )
		{
			types::fcolor color_src( image_source.At( x, y ) );
			types::fcolor color_dest( dest.At( x, y ) );

			color_dest.r *= color_src.r ;
			color_dest.g *= color_src.g;
			color_dest.b *= color_src.b ;

			color_dest.a = 1.f;

			dest.At( x, y ) = color_dest.Get32();
		}
	}
}

bool IsGoodMarker( const Marker& marker )
{
	for( int i = 0; i < marker.corners.length; ++i )
	{
		for( int j = i + 1; j < marker.corners.length; ++j )
		{
			if( marker.corners[i] == marker.corners[j] )
				return false;
		}
	}

	float l1 = ( marker.corners[C_TOP_LEFT] - marker.corners[C_TOP_RIGHT] ).Length();
	float l2 = ( marker.corners[C_BOTTOM_LEFT] - marker.corners[C_BOTTOM_RIGHT] ).Length();
	
	float diff = abs( l1 - l2 );
	float diff_p = diff / l1;
	if( diff_p > 0.2f ) return false;

	l1 = ( marker.corners[C_TOP_LEFT] - marker.corners[C_BOTTOM_LEFT] ).Length();
	l2 = ( marker.corners[C_TOP_RIGHT] - marker.corners[C_BOTTOM_RIGHT] ).Length();
	
	diff = abs( l1 - l2 );
	diff_p = diff / l1;
	if( diff_p > 0.2f ) return false;

	return true;
}

//============================================================================================


void DoCard( const std::string& input_file, const std::string& card_image, const std::string& output_filename )
{
	ParseQR_LoadImage( input_file );	

	bool is_good = false;

	for( int iteration = 0; iteration < 12 && is_good == false; ++iteration )
	{
		mLines.clear();
		mMarkers.clear();
		mDebugMarkers.clear();

		ParseQR( input_file, iteration );

	#if 1
		// check for once that intersect with other
		std::vector< LinePair > new_lines;
		
		using namespace ceng::math;

		// keep x longest
		{
			float shortest = 0;
			int how_many = 48;
			for( std::size_t i = 0; i < mLines.size(); ++i )
			{
				float length = Distance( mLines[i].a.a, mLines[i].a.b );
				if( length > shortest )
				{
					new_lines.push_back( mLines[i] );
					
					if( (int)new_lines.size() > how_many )
					{
						RemoveShortest( new_lines );

						shortest = FLT_MAX;
						int shortest_i = -1;
						for( int i = 0; i < (int)new_lines.size(); ++i )
						{
							float l = Distance( new_lines[i].a.a, new_lines[i].a.b );
							if( l < shortest )
							{
								shortest_i = i;
								shortest = l;
							}
						}
					}
				}
			}

			mLines = new_lines;
			new_lines.clear();
		}

		// keep only the ones that intersect
		{
			std::vector< bool > lines_used( mLines.size() );
			for( std::size_t i = 0; i < lines_used.size(); ++i )
				lines_used[i] = false;

			for( std::size_t i = 0; i < mLines.size(); ++i )
			{
				LinePair li = mLines[i];

				bool found_a_pair = false;
				Vec2 result;
				for( std::size_t j = i + 1; j < mLines.size(); ++j )
				{
					LinePair lj = mLines[j];

					if( LineIntersection( li.a.a, li.b.a, lj.a.a, lj.b.a, result ) &&
						LineIntersection( li.a.a, li.b.a, lj.a.b, lj.b.b, result ) &&
						LineIntersection( li.a.b, li.b.b, lj.a.a, lj.b.a, result ) &&
						LineIntersection( li.a.b, li.b.b, lj.a.b, lj.b.b, result ) )
					{
						// these are fine...
						// AddAsMarker( li, lj );
						found_a_pair = true;
						if( lines_used[j] == false )
						{
							lines_used[j] = true;
							new_lines.push_back( lj );
						}
					}
				}

				if( found_a_pair && lines_used[i] == false )
				{
					lines_used[i] = true;
					new_lines.push_back( li );
				}

			}
		}

		mLines = new_lines;
		new_lines.clear();

		// look for 8 same lengthish
		{
			const int how_many = 4 * 2 - 1;
			std::vector< float > line_l( mLines.size() );
			std::vector< float > sorted_line_l( mLines.size() );

			std::vector< float > lowest_ls;
			float lowest_value = FLT_MAX;

			for( std::size_t i = 0; i < mLines.size(); ++i )
			{
				float length = Distance( mLines[i].a.a, mLines[i].b.a );
				if( length == 0 ) 
					continue;
				
				for( std::size_t j = 0; j < mLines.size(); ++j )
				{
					float l = Distance( mLines[j].a.a, mLines[j].b.a );
					float lp = l / length;
					line_l[j] = abs( 1.f - lp );
				}

				// check which are the 8 sameish lengths
				sorted_line_l = line_l;
				std::sort( sorted_line_l.begin(), sorted_line_l.end() );
				
				// check the 8th most different length
				if( sorted_line_l.size() > how_many && sorted_line_l[ how_many ] < lowest_value )
				{
					lowest_value = sorted_line_l[ how_many ];
					lowest_ls = line_l;
				}
			}

			// check if the lowest value is low enough...
			// now we should have the lowest_value 
			for( std::size_t i = 0; i < lowest_ls.size(); ++i )
			{
				if( lowest_ls[i] <= lowest_value + 0.21f )
				{
					new_lines.push_back( mLines[i] );
				}
			}

			
			mLines = new_lines;
			new_lines.clear();
		}

		// look for the markers
		for( std::size_t i = 0; i < mLines.size(); ++i )
		{
			LinePair li = mLines[i];

			bool found_a_pair = false;
			Vec2 result;
			for( std::size_t j = i + 1; j < mLines.size(); ++j )
			{
				LinePair lj = mLines[j];

				if( LineIntersection( li.a.a, li.b.a, lj.a.a, lj.b.a, result ) &&
					LineIntersection( li.a.a, li.b.a, lj.a.b, lj.b.b, result ) &&
					LineIntersection( li.a.b, li.b.b, lj.a.a, lj.b.a, result ) &&
					LineIntersection( li.a.b, li.b.b, lj.a.b, lj.b.b, result ) )
				{
					// these are fine...
					AddAsMarker( li, lj );
					found_a_pair = true;
					new_lines.push_back( lj );
				}
			}

			if( found_a_pair ) 
				new_lines.push_back( li );

		}

		mLines = new_lines;
		new_lines.clear();

		// remove the wonky markers
		for( std::size_t i = 0; i < mMarkers.size(); )
		{
			if( IsGoodMarker( mMarkers[i] ) )
			{
				++i;
			}
			else
			{
				mMarkers[i] = mMarkers[ mMarkers.size() - 1 ];
				mMarkers.pop_back();
			}
		}

		// now let's look if we have 4 corners
		if( mMarkers.size() >= 4 )
		{
			std::cout << "markers: " << mMarkers.size() << std::endl;
			mDebugMarkers = mMarkers;


			types::vector2 aabb_min( FLT_MAX, FLT_MAX );
			types::vector2 aabb_max( FLT_MIN, FLT_MIN );

			for( std::size_t i = 0; i < mMarkers.size(); ++i )	
			{
				mMarkers[i].UpdateAABB();
				if( mMarkers[i].center.x < aabb_min.x ) aabb_min.x = mMarkers[i].center.x;
				if( mMarkers[i].center.y < aabb_min.y ) aabb_min.y = mMarkers[i].center.y;
				if( mMarkers[i].center.x > aabb_max.x ) aabb_max.x = mMarkers[i].center.x;
				if( mMarkers[i].center.y > aabb_max.y ) aabb_max.y = mMarkers[i].center.y;
			}

			std::vector< Marker > markers( 4 );
			markers[C_TOP_LEFT] = FindClosestTo( mMarkers,	types::vector2( aabb_min.x, aabb_min.y ) );
			markers[C_TOP_RIGHT] = FindClosestTo( mMarkers, types::vector2( aabb_max.x, aabb_min.y ) );
			markers[C_BOTTOM_RIGHT] = FindClosestTo( mMarkers, types::vector2( aabb_max.x, aabb_max.y ) );
			markers[C_BOTTOM_LEFT] = FindClosestTo( mMarkers, types::vector2( aabb_min.x, aabb_max.y ) );
			
			mMarkers = markers;
		}

		if( mMarkers.size() == 4 )
		{
			bool overlaps = false;
			
			for( std::size_t i = 0; i < mMarkers.size(); ++i )
			{
				for( std::size_t j = i + 1; j < mMarkers.size(); ++j )
				{
					if( ceng::math::TestAABBAABBIntersection( 
						mMarkers[i].aabb_min, mMarkers[i].aabb_max,
						mMarkers[j].aabb_min, mMarkers[j].aabb_max ) )
					{
						overlaps = true;
					}
				}
			}

			if( overlaps == false )
			{
				is_good = true;
				break;
			}
			else
			{
				std::cout << "found overlaps" << std::endl;
			}
		}
	}

	// -------------- now if we have the marks ------

	/*ceng::CArray2D< unsigned int > image_texture;
	image_texture = mImage;*/
	// LoadImageTo( "test/test_012.png", image_texture );	

	ceng::CArray2D< unsigned int > temp_texture(mImage.GetWidth(), mImage.GetHeight());
	// ceng::CArray2D< unsigned int > temp_texture2(mImage.GetWidth(), mImage.GetHeight());


	// Fill the markers
	for( int y = 0; y < temp_texture.GetHeight(); ++y )
	{
		for( int x = 0; x < temp_texture.GetWidth(); ++x )
		{
			if(  temp_texture.At( x, y ) == 0 )
			{
				types::fcolor color;
				color.Set32( mImage.At( x, y ) );
				color.SetA( 0 );
				temp_texture.At( x, y ) = color.Get32();
			}
		}
	}

	for( std::size_t i = 0; i < mMarkers.size(); ++i )
	{
		FillTheBlanks( temp_texture, mMarkers[i] );
	}

	for( std::size_t i = 0; i < mMarkers.size(); ++i )
	{
		Blur( 1.0f, temp_texture, mMarkers[i].aabb_min - types::vector2( 3, 3 ), mMarkers[i].aabb_max + types::vector2( 3, 3 ) );
		BlitTo( temp_texture, mImage, 
			(int)mMarkers[i].aabb_min.x - 7,
			(int)mMarkers[i].aabb_min.y - 7,
			(int)mMarkers[i].aabb_max.x + 7,
			(int)mMarkers[i].aabb_max.y + 7 );
	}

	// 

	// test image
	if( mMarkers.size() >= 4 )
	{
		// "c8.png"
		ceng::CArray2D< unsigned int > card_texture;
		LoadImageTo( card_image, card_texture );

		temp_texture.SetEverythingTo( 0xFFFFFFFF );

		float rate = 2.f / 10.f;
		ExpandCorner( mMarkers[C_TOP_LEFT], C_TOP_LEFT, rate );
		ExpandCorner( mMarkers[C_TOP_RIGHT], C_TOP_RIGHT, rate );
		ExpandCorner( mMarkers[C_BOTTOM_LEFT], C_BOTTOM_LEFT, rate );
		ExpandCorner( mMarkers[C_BOTTOM_RIGHT], C_BOTTOM_RIGHT, rate );

		types::vector2 aabb_min( FLT_MAX, FLT_MAX );
		types::vector2 aabb_max( FLT_MIN, FLT_MIN );
		for( std::size_t i = 0; i < mMarkers.size(); ++i )
		{
			aabb_min.x = ceng::math::Min( aabb_min.x, mMarkers[i].aabb_min.x );
			aabb_min.y = ceng::math::Min( aabb_min.y, mMarkers[i].aabb_min.y );
			aabb_max.x = ceng::math::Max( aabb_max.x, mMarkers[i].aabb_max.x );
			aabb_max.y = ceng::math::Max( aabb_max.y, mMarkers[i].aabb_max.y );
		}

		types::vector2 center_p(0,0);
		for( int i = 0; i < 4; ++i )
			center_p += mMarkers[i].corners[i];

		center_p = (1.f / 4.f) * center_p;

		std::vector< types::vector2 > triangle(3);
		std::vector< types::vector2 > text_coords(3);

		// top
		triangle[0] = ( mMarkers[C_TOP_LEFT].corners[C_TOP_LEFT] );
		triangle[1] = ( mMarkers[C_TOP_RIGHT].corners[C_TOP_RIGHT] );
		triangle[2] = center_p;

		text_coords[0].Set( 0, 0 );
		text_coords[1].Set( (float)card_texture.GetWidth(), 0 );
		text_coords[2].Set( 0.5f * (float)card_texture.GetWidth(), 0.5f * (float)card_texture.GetHeight() );

		DrawTriangle( card_texture, text_coords, triangle, temp_texture );

		// right
		triangle[0] = ( mMarkers[C_TOP_RIGHT].corners[C_TOP_RIGHT] );
		triangle[1] = ( mMarkers[C_BOTTOM_RIGHT].corners[C_BOTTOM_RIGHT] );
		triangle[2] = center_p;

		text_coords[0].Set( (float)card_texture.GetWidth(), 0 );
		text_coords[1].Set( (float)card_texture.GetWidth(), (float)card_texture.GetHeight() );
		text_coords[2].Set( 0.5f * (float)card_texture.GetWidth(), 0.5f * (float)card_texture.GetHeight() );

		DrawTriangle( card_texture, text_coords, triangle, temp_texture );

		// bottom
		triangle[0] = center_p;
		triangle[1] = ( mMarkers[C_BOTTOM_RIGHT].corners[C_BOTTOM_RIGHT] );
		triangle[2] = ( mMarkers[C_BOTTOM_LEFT].corners[C_BOTTOM_LEFT] );

		text_coords[0].Set( 0.5f * (float)card_texture.GetWidth(), 0.5f * (float)card_texture.GetHeight() );
		text_coords[1].Set( (float)card_texture.GetWidth(), (float)card_texture.GetHeight() );
		text_coords[2].Set( 0, (float)card_texture.GetHeight() );

		DrawTriangle( card_texture, text_coords, triangle, temp_texture );

		// left
		triangle[0] = center_p;
		triangle[2] = ( mMarkers[C_BOTTOM_LEFT].corners[C_BOTTOM_LEFT] );
		triangle[1] = ( mMarkers[C_TOP_LEFT].corners[C_TOP_LEFT] );

		text_coords[0].Set( 0.5f * (float)card_texture.GetWidth(), 0.5f * (float)card_texture.GetHeight() );
		text_coords[1].Set( 0, (float)card_texture.GetHeight() );
		text_coords[1].Set( 0, 0 );

		DrawTriangle( card_texture, text_coords, triangle, temp_texture );


#if 0
		// old, wrong perpective
		std::vector< types::vector2 > triangle(3);

		triangle[0] = ( mMarkers[C_TOP_LEFT].corners[C_TOP_LEFT] );
		triangle[1] = ( mMarkers[C_TOP_RIGHT].corners[C_TOP_RIGHT] );
		triangle[2] = ( mMarkers[C_BOTTOM_RIGHT].corners[C_BOTTOM_RIGHT] );

		std::vector< types::vector2 > text_coords(3);
		text_coords[0].Set( 0, 0 );
		text_coords[1].Set( (float)card_texture.GetWidth(), 0 );
		text_coords[2].Set( (float)card_texture.GetWidth(), (float)card_texture.GetHeight() );

		DrawTriangle( card_texture, text_coords, triangle, temp_texture );


		triangle[0] = ( mMarkers[C_TOP_LEFT].corners[C_TOP_LEFT] );
		triangle[1] = ( mMarkers[C_BOTTOM_RIGHT].corners[C_BOTTOM_RIGHT] );
		triangle[2] = ( mMarkers[C_BOTTOM_LEFT].corners[C_BOTTOM_LEFT] );

		text_coords[0].Set( 0, 0 );
		text_coords[1].Set( (float)card_texture.GetWidth(), (float)card_texture.GetHeight() );
		text_coords[2].Set( 0, (float)card_texture.GetHeight() );

		DrawTriangle( card_texture, text_coords, triangle, temp_texture );
#endif
		// Blur( 1.f, temp_texture2, aabb_min - types::vector2( 3, 3 ), aabb_max + types::vector2( 3, 3 ) );
		GaussBlur( 1.3f, temp_texture,
			(int)aabb_min.x - 10, 
			(int)aabb_min.y - 10,
			(int)aabb_max.x + 10,
			(int)aabb_max.y + 10 );

		BlitMultiply( temp_texture, mImage, 
			(int)aabb_min.x - 7, 
			(int)aabb_min.y - 7,
			(int)aabb_max.x + 7,
			(int)aabb_max.y + 7 ); 
	} 

#endif
	
	// height
	int new_height = (int)( ( 800.f / (float)mImage.GetWidth() ) * (float)mImage.GetHeight() + 0.5f );
	ResizeImage( mImage, 800, new_height );

	SaveImageTo( output_filename, mImage );
}



#ifdef COMMAND_LINE

int main(int argc, char *argv[])
{

	if( argc < 4 )
	{
		std::cout << "usage: " << std::endl;
		std::cout << "parse.exe source_image.jpg card_image.png output_image.jpg" << std::endl;
		return 0;
	}

	srand(time(NULL));

	std::string source_image = argv[1];
	std::string card_image = argv[2];
	std::string output_image = argv[3];

	DoCard( source_image, card_image, output_image );

	// std::cout << "hello world" << std::endl;
	// ParseQR( "test/test_006.png" );
	return 0;
}

#endif