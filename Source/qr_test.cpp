#include "qr_test.h"

// #define IMPL_DEV_BUILD

#include <game_utils/drawlines/drawlines.h>

#include <utils/staticarray/cstaticarray.h>
#include <utils/math/cstatisticshelper.h>
#include <utils/vector_utils/vector_utils.h>
#include <utils/color/ccolor.h>
#include <utils/math/point_inside.h>

#include "gameplay_utils/game_mouse.h"
#include "misc_utils/debug_layer.h"
#include "misc_utils/simple_profiler.h"

#include "qr_parser.cpp"


// ----------------------------------------------------------------------------

bool config_display_wireframe = false;


QRTest::QRTest()
{
}


void QRTest::Exit()
{
	std::cout << "QRTest::EXit() called" << std::endl;
	mDebugLayer.reset( NULL );
}

//-----------------------------------------------------------------------------
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

std::vector< Marker > mMarkers;

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
	ExtendLine( top_start, top_end, 1024.f );

	vector2 bottom_start( l0.b.a );
	vector2 bottom_end( l0.b.b );
	ExtendLine( bottom_start, bottom_end, 1024.f );

	vector2 left_start( l1.a.a );
	vector2 left_end( l1.a.b );
	ExtendLine( left_start, left_end, 1024.f );

	vector2 right_start( l1.b.a );
	vector2 right_end( l1.b.b );
	ExtendLine( right_start, right_end, 1024.f );

	types::vector2 top_center = Lerp( top_start, top_end, 0.5f );
	types::vector2 bottom_center = Lerp( bottom_start, bottom_end, 0.5f );
	types::vector2 left_center = Lerp( left_start, left_end, 0.5f );
	types::vector2 right_center = Lerp( right_start, right_end, 0.5f );

	// look for the intersection
	if( LineIntersection( top_start, top_end, left_start, left_end, result.corners[C_TOP_LEFT] ) == false )
	{
		std::cout << "Error couldn't find top left for marker..." << std::endl;
		return false;
	}

	if( LineIntersection( top_start, top_end, right_start, right_end, result.corners[C_TOP_RIGHT] ) == false )
	{
		std::cout << "Error couldn't find top left for marker..." << std::endl;
		return false;
	}

	if( LineIntersection( bottom_start, bottom_end, right_start, right_end, result.corners[C_BOTTOM_RIGHT] ) == false )
	{
		std::cout << "Error couldn't find top left for marker..." << std::endl;
		return false;
	}

	if( LineIntersection( bottom_start, bottom_end, left_start, left_end, result.corners[C_BOTTOM_LEFT] ) == false )
	{
		std::cout << "Error couldn't find top left for marker..." << std::endl;
		return false;
	}
	
	if( LineIntersection( top_center, bottom_center, left_center, right_center, result.center ) == false )
	{
		std::cout << "Error couldn't find center for marker" << std::endl;
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
							c = image_data[ pos_y * image.GetWidth() + pos_x ];
						}
						else
						{
							int pos_x = int( ( ceng::Randomf( 0.f, 1.f ) > ix ) ? aabb_min.x - 2 : aabb_max.x + 2 );
							int pos_y = y + ceng::Random( -5, 5 );
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


//============================================================================================

void DoCard( const std::string& input_file, const std::string& card_image, const std::string& output_filename )
{
	ParseQR_LoadImage( input_file );	

	bool is_good = false;

	for( int iteration = 0; iteration < 12 && is_good == false; ++iteration )
	{
		ParseQR( input_file, iteration );

	#if 1
		// check for once that intersect with other
		std::vector< LinePair > new_lines;
		
		using namespace ceng::math;

		// keep x longest
		{
			float shortest = 0;
			int how_many = 16;
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
				if( lowest_ls[i] <= lowest_value + 0.1f )
				{
					new_lines.push_back( mLines[i] );
				}
			}

			
			mLines = new_lines;
			new_lines.clear();
		}
		
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

		// now let's look if we have 4 corners
		if( mMarkers.empty() == false )
		{
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
	
	ResizeImage( mImage, 800, 600 );

	SaveImageTo( output_filename, mImage );
}


void QRTest::Init()
{
	DefaultApplication::Init();
	Poro()->GetGraphics()->SetFillColor( poro::GetFColor( 248.f / 255.f, 245.f / 255.f, 236.f / 255.f, 1.f ) );

	mSpriteContainer = new as::Sprite;
	mDebugLayer.reset( new DebugLayer );

	DoCard( "test/t1.jpg", "cards/h12.png", "test_image.jpg" );

	// --- graphics for poro ---
	as::Sprite* sprite = as::LoadSprite( "test_image.jpg" );
	sprite->SetScale( 1024.f / sprite->GetSize().x, 768.f / sprite->GetSize().y );
	mSpriteContainer->addChild( sprite );
	

}

// ----------------------------------------------------------------------------

void QRTest::Update( float dt )
{
	// UpdateGTweens( dt );

	if( mDebugLayer.get() ) 
		mDebugLayer->Update( dt );

	if( mSpriteContainer )
		mSpriteContainer->Update( dt );

	GameMouse::GetSingletonPtr()->OnFrameEnd();
}

// ----------------------------------------------------------------------------

void QRTest::Draw( poro::IGraphics* graphics )
{ 
	SetLineWidth( 3.f );

	SPROFILE( "QRTest::Draw()" );

	if( mDebugLayer.get() ) 
		mDebugLayer->Draw( graphics );


	SetLineWidth( 2.f );


	if( mSpriteContainer )
		as::DrawSprite( mSpriteContainer, graphics );
	

#if 0
	for( std::size_t i = 0; i < mLines.size(); ++i )
	{
		DrawLine( graphics, 
			types::vector2( mLines[i].a.a ),  
			types::vector2( mLines[i].a.b ),
			poro::GetFColor( 1, 0, 0, 1.0f ),
			NULL );

		DrawLine( graphics, 
			types::vector2( mLines[i].b.a ),  
			types::vector2( mLines[i].b.b ),
			poro::GetFColor( 1, 0, 0, 1.0f ),
			NULL );

		// maybe draw the paraller lines
		DrawLine( graphics, 
			types::vector2( mLines[i].a.a ),  
			types::vector2( mLines[i].b.a ),
			poro::GetFColor( 0, 1, 0, 0.5f ),
			NULL );

		DrawLine( graphics, 
			types::vector2( mLines[i].a.b ),  
			types::vector2( mLines[i].b.b ),
			poro::GetFColor( 0, 1, 0, 0.5f ),
			NULL );
		
	}
#endif

	// 
	if( config_display_wireframe && mMarkers.empty() == false )
	{
		types::vector2 center_p(0,0);
		for( std::size_t k = 0; k < mMarkers.size(); ++k )	
		{
			center_p += mMarkers[k].center;
			int j, i;
			for (j=4-1,i=0; i < 4; j=i++) 
			{
				DrawLine( graphics, 
					( mMarkers[k].corners[i] ),  
					( mMarkers[k].corners[j] ),
					poro::GetFColor( 0, 1, 0, 0.5f ),
					NULL );
			}
		}


		center_p.x /= (float)mMarkers.size();
		center_p.y /= (float)mMarkers.size();

		DrawCircle( graphics, center_p, 7.f, poro::GetFColor( 0, 0, 1, 1 ) );

		// -- draw the quad
		DrawLine( graphics, 
				( mMarkers[C_TOP_LEFT].corners[C_TOP_LEFT] ),  
				( mMarkers[C_TOP_RIGHT].corners[C_TOP_RIGHT] ),
				poro::GetFColor( 1, 0, 0, 0.5f ),
				NULL );

		DrawLine( graphics, 
			( mMarkers[C_TOP_RIGHT].corners[C_TOP_RIGHT] ),  
			( mMarkers[C_BOTTOM_RIGHT].corners[C_BOTTOM_RIGHT] ),
			poro::GetFColor( 1, 0, 0, 0.5f ),
			NULL );

		DrawLine( graphics, 
			( mMarkers[C_BOTTOM_RIGHT].corners[C_BOTTOM_RIGHT] ),  
			( mMarkers[C_BOTTOM_LEFT].corners[C_BOTTOM_LEFT] ),
			poro::GetFColor( 1, 0, 0, 0.5f ),
			NULL );

		DrawLine( graphics, 
			( mMarkers[C_BOTTOM_LEFT].corners[C_BOTTOM_LEFT] ),  
			( mMarkers[C_TOP_LEFT].corners[C_TOP_LEFT] ),
			poro::GetFColor( 1, 0, 0, 0.5f ),
			NULL );
	}
}

// ----------------------------------------------------------------------------

void QRTest::MouseMove(const poro::types::vec2& p)
{
	if (mDebugLayer.get())
		mDebugLayer->mMousePositionInWorld = types::vector2(p);
}

void QRTest::MouseButtonDown(const poro::types::vec2& p, int button)
{
}

void QRTest::MouseButtonUp(const poro::types::vec2& pos, int button)
{
}

//=============================================================================

void QRTest::OnKeyDown( int key, poro::types::charset unicode )
{
	if( key == 27 ) 
		Poro()->Exit();

	// p
	if( key == 112 && mSpriteContainer )
	{
		as::Sprite* overlay = mSpriteContainer->GetChildByName( "overlay" );
		if( overlay )
			overlay->SetVisibility( !overlay->GetVisibility() );

	}

	// w
	if( key == 119 )
	{
		config_display_wireframe = !config_display_wireframe;
	}

}

//-----------------------------------------------------------------------------


void QRTest::OnKeyUp( int key, poro::types::charset unicode )
{
}

//=============================================================================

