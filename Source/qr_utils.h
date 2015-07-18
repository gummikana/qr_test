#ifndef INC_QR_UTILS_H
#define INC_QR_UTILS_H

#include <time.h>
#include <cstdlib>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#define cassert assert


namespace ceng {
namespace math {


template< typename Type >
inline const Type& Min( const Type& c1, const Type& c2 ) {
	return ( c1 < c2 )?c1:c2;
}

template< typename Type >
inline const Type& Max( const Type& c1, const Type& c2 ) {
	return ( c1 < c2 )?c2:c1;
}

template< typename Type >
inline Type Clamp( const Type& a, const Type& low, const Type& high) {
	return ceng::math::Max( low, ceng::math::Min( a, high ) );
}


template< typename Type >
inline Type Lerp( const Type& a, const Type& b, float t ) {
	return (Type)( a + t * (b - a) );
}

template< class T >
bool LineIntersection( const T& startA, const T& endA,
	const T& startB, const T& endB, T& result )
{
	//TODO: reuse mathutil.intersect
	float d = (float)( (endB.y - startB.y) * (endA.x - startA.x) - (endB.x - startB.x) * (endA.y - startA.y) );
	
	if ( d == 0 ) // parallel lines
		return false;
	
	float uA = (float)( (endB.x - startB.x) * (startA.y - startB.y) - (endB.y - startB.y) * (startA.x - startB.x) );
	uA /= d;
	float uB = (float)( (endA.x - startA.x) * (startA.y - startB.y) - (endA.y - startA.y) * (startA.x - startB.x) );
	uB /= d;
	
	if ( uA < 0 || uA > 1 || uB < 0 || uB > 1 ) 
		return false; // intersection point isn't between the start and endpoints


	result.x = (startA.x + uA * (endA.x - startA.x) );
	result.y = (startA.y + uA * (endA.y - startA.y) );

	// result = position;
	return true;
}

// -- vector2 ---

template< class Type >
class CVector2
{
public:
	typedef Type unit_type;
	//=========================================================================

	CVector2() : x( Type() ), y( Type() ) { }

	CVector2( Type x, Type y) : x( x ), y( y ) { }

	template< class T >
	explicit CVector2( const T& other ) : x( (Type)( other.x ) ), y( (Type)(other.y) ) { }


	void Set( Type x_, Type y_) { x = x_; y = y_; }

	//=========================================================================

	bool operator== ( const CVector2< Type >& other ) const { return ( this->x == other.x && this->y == other.y ); }
	bool operator!= ( const CVector2< Type >& other ) const { return !operator==( other ); }

	//=========================================================================

	CVector2< Type > operator -() const { return CVector2< Type >( -x, -y ); }
	
	CVector2< Type >& operator += ( const CVector2< Type >& v )
	{
		x += v.x; y += v.y;

		return *this;
	}
	
	CVector2< Type >& operator -= ( const CVector2< Type >& v )
	{
		x -= v.x; y -= v.y;

		return *this;
	}

	
	CVector2< Type >& operator *= ( float a )
	{
		x = (Type)(x * a); y = (Type)(y * a);

		return *this;
	}
	
	CVector2< Type >& operator /= ( float a )
	{
		x = (Type)(x / a); y = (Type)(y / a);

		return *this;
	}

	//-------------------------------------------------------------------------
	
	CVector2< Type > operator + ( const CVector2< Type >& other ) const
	{
		return CVector2< Type >( this->x + other.x, this->y + other.y );
	}

	CVector2< Type > operator - ( const CVector2< Type >& other ) const
	{
		return CVector2< Type >( this->x - other.x, this->y - other.y );
	}

	CVector2< Type > operator * ( float t ) const
	{
		return CVector2< Type >( (Type)( this->x * t ), (Type)( this->y * t ) );
	}

	CVector2< Type > operator / ( float t ) const
	{
		return CVector2< Type >( (Type)( this->x / t ), (Type)( this->y / t ) );
	}
	
	//=========================================================================

	Type Dot( const CVector2< Type >& a, const CVector2< Type >& b ) const
	{
		return a.x * b.x + a.y * b.y;
	}

	Type Cross( const CVector2< Type >& a, const CVector2< Type >& b ) const
	{
		return a.x * b.y - a.y * b.x;
	}

	CVector2< Type > Cross( const CVector2< Type >& a, const Type& s ) const
	{
		return CVector2< Type >( s * a.y, -s * a.x );
	}

	CVector2< Type > Cross( const Type& s, const CVector2< Type >& a ) const
	{
		return CVector2< Type >( -s * a.y, s * a.x );
	}

	//=========================================================================

	Type LengthSquared() const
	{
		return ( x * x + y * y );
	}

	float Length() const
	{
		return sqrtf( (float)LengthSquared() );
	}

	CVector2< Type > Normalize() const
	{
		float d = Length();
		
		if( d > 0 )
			return CVector2< Type >( Type( x / d ), Type( y / d ) );
		else 
			return CVector2< Type >( 0, 0 );
	}

	//=========================================================================
	// ripped from:
	// http://forums.indiegamer.com/showthread.php?t=10459

	Type Angle() const 
	{
		CVector2< Type > normal = Normalize();
		Type angle = atan2( normal.y, normal.x );

		return angle;
	}

	Type Angle( const CVector2< Type >& x ) const
	{
		Type dot = Dot( *this, x );
		Type cross = Cross( *this, x );
		
		// angle between segments
		Type angle = (Type) atan2( cross, dot );

		return angle;
	}

	CVector2< Type >& Rotate( Type angle_rad ) 
	{
		Type tx = x;
		x = (Type)x * (Type)cos(angle_rad) - y * (Type)sin(angle_rad);
		y = (Type)tx * (Type)sin(angle_rad) + y * (Type)cos(angle_rad);

		return *this;
	}

	CVector2< Type >& Rotate( const CVector2< Type >& centre, Type angle_rad )
	{
		CVector2< Type > D = *this - centre;
		D.Rotate( angle_rad );

		// *this = xCentre + D;
		D += centre;
		Set( D.x, D.y );

		return *this;
	}

	//=========================================================================

	Type x;
	Type y;
};

template< class Type >
CVector2< Type > operator * ( Type s, const CVector2< Type >& v) {
	return CVector2< Type >(s * v.x, s * v.y);
}

// ---------------------------------------------------------------------------------------

template< class PType >
bool TestAABBAABBIntersection( const CVector2< PType >& a_min, const CVector2< PType >& a_max,
							  const CVector2< PType >& b_min, const CVector2< PType >& b_max )
{
    if (a_max.x < b_min.x || 
        a_max.y < b_min.y || 
        a_min.x > b_max.x || 
        a_min.y > b_max.y) 
    {
        return false;
    }
    return true;
}


} // end  of namespace math
} // end of namespace ceng

// ---------- types ----------

namespace ceng {


template< class T, int N >
class CStaticArray
{
public:
	CStaticArray() :
		data(),
		length( N )
	{
		for( int i = 0; i< N; ++i )
		{
			data[ i ] = T();
		}
	}

	CStaticArray( const CStaticArray< T, N >& other ) :
		data(),
		length( N )
	{	for( int i = 0; i < N; ++i )
		{
			data[ i ] = other.data[ i ];
		}
	}

	~CStaticArray()
	{
		// delete [] data;
	}

	CStaticArray& operator=( const CStaticArray< T, N >& other )
	{
		for( int i = 0; i < N; ++i )
		{
			data[ i ] = other.data[ i ];
		}

		return *this;
	}


	T& operator[] ( int i )
	{
		cassert( i >= 0 && i < N );
		return data[ i ];
	}

	const T& operator[]( int i ) const
	{
		cassert( i >= 0 && i < N );
		return data[ i ];
	}

	const int length;

private:

	T data[ N ];

};

// --- safe array


template< typename Type, typename SizeType = int >
class CSafeArray
{
public:
	CSafeArray() : data( 0 ), _size( SizeType() ) { }
	CSafeArray( SizeType size ) : 
		data( new Type[ size ] ),
		_size( size )
	{
		
		/*for( SizeType i = 0; i < Size(); ++i ) 
			data[ i ] = Type();*/
		
		memset( data, (int)Type(), Size() * sizeof( Type ) );
	}
	
	CSafeArray( const CSafeArray& other ) :
		data( 0 ), 
		_size( SizeType() )
	{
		operator=(other);
	}

	~CSafeArray()
	{
		Clear();
	}

	const CSafeArray& operator=( const CSafeArray& other )
	{
		if( other._size != _size )
		{
			Clear();
			data =  new Type[ other._size ];
			_size = other._size;
		}

		/*
		for( SizeType i = 0; i < Size(); ++i ) 
			data[ i ] = other.Rand( i );
		*/
		memcpy( data, other.data, Size() * sizeof( Type ) );
		// fast_memcpy( data, other.data, Size() * sizeof( Type ) );
		// X_aligned_memcpy_sse2( data, other.data, Size() * sizeof( Type ) );

		return *this;
	}

	inline Type& operator[]( SizeType i )
	{
		cassert( !( i < 0 || i >= _size ) );
		return data[ i ];
	}

	inline const Type& operator[]( SizeType i ) const 
	{
		cassert( !( i < 0 || i >= _size ) );
		return data[ i ];
	}

	void Clear()
	{
		delete [] data;
		data = 0;
		_size = 0;
	}

	void clear() { Clear(); }

	SizeType Size() const { return _size; }
	SizeType size() const { return _size; }

	bool Empty() const { return _size == 0; }
	bool empty() const { return Empty(); }

	void Resize( SizeType s ) 
	{
		if( _size != s )
		{
			Clear();
			data =  new Type[ s ];
			_size = s;
			
			
			/*for( SizeType i = 0; i < Size(); ++i ) 
				data[ i ] = Type();*/
			memset( data, (int)Type(), Size() * sizeof( Type ) );
		}
	}

	void resize( SizeType s ) { Resize( s ); }

	inline const Type& At( SizeType i ) const
	{
		if( i < 0 || i >= _size )
			return Type();

		return data[ i ];
	}

	inline const Type& Rand( SizeType i ) const
	{
		cassert( !( i < 0 || i >= _size ) );

		return data[ i ];
	}

	inline Type& Rand( SizeType i )
	{
		cassert( !( i < 0 || i >= _size ) );

		return data[ i ];
	}

	Type* data;
private:
	SizeType _size;
};

// --- carray2d ---


template < class _Ty, class _A = std::allocator<_Ty> >
class CArray2D
{
public:

	typedef typename _A::reference		reference;
	typedef typename _A::const_reference const_reference;

	class CArray2DHelper
	{
	public:
		CArray2DHelper( CArray2D& array ) : myArray( array ) { }
		~CArray2DHelper() { }

		reference operator [] ( int _y ) { return myArray.Rand( myX, _y ); }
		const_reference operator [] ( int _y ) const { return myArray.Rand( myX, _y ); }

		void SetX( int _x ) const { myX = _x; }

	private:
		mutable int			myX;
		CArray2D&	myArray;

	};

	CArray2D() :
	  myWidth( 0 ),
	  myHeight( 0 ),
	  mySize( 0 ),
	  myArraysLittleHelper( *this ),
	  myNullReference( _Ty() )
	{

	}

	CArray2D( int _width, int _height ) :
	  myWidth( _width ),
	  myHeight( _height ),
	  mySize( 0 ),
	  myArraysLittleHelper( *this ),
	  myNullReference( _Ty() )
	{

		Allocate();
	}

	CArray2D( const CArray2D< _Ty, _A >& other ) :
		myWidth( other.myWidth ),
		myHeight( other.myHeight ),
		mySize( other.mySize ),
		myArraysLittleHelper( *this ),
		myDataArray( other.myDataArray ),
	  myNullReference( _Ty() )
	{

	}

	CArray2DHelper& operator[] ( int _x ) { myArraysLittleHelper.SetX( _x ); return myArraysLittleHelper; }
	const CArray2DHelper& operator[] ( int _x ) const { myArraysLittleHelper.SetX( _x ); return myArraysLittleHelper; }

	const CArray2D< _Ty, _A >& operator=( const CArray2D< _Ty, _A >& other )
	{
		myWidth = other.myWidth;
		myHeight = other.myHeight;
		mySize = other.mySize;
		myDataArray = other.myDataArray;

		return *this;
	}

	inline int GetWidth() const
	{
		return myWidth;
	}

	inline int GetHeight() const
	{
		return myHeight;
	}

	inline bool IsValid( int x, int y ) const 
	{
		if( x < 0 ) return false;
		if( y < 0 ) return false;
		if( x >= myWidth ) return false;
		if( y >= myHeight ) return false;
		return true;
	}

	void Resize( int _width ) { SetWidth( _width ); }
	void Resize( int _width, int _height ) { SetWidthAndHeight( _width, _height ); }

	void SetWidth(  int _width  ) { myWidth  = _width; Allocate(); }
	void SetHeight( int _height ) { myHeight = _height; Allocate(); }

	void SetWidthAndHeight( int _width, int _height ) { myWidth = _width; myHeight = _height; Allocate(); }

	void SetEverythingTo( const _Ty& _who )
	{
		int i;
		for ( i = 0; i < mySize; i++ )
			myDataArray[ i ] = _who;
	}


	inline reference At( int _x, int _y )
	{
#ifdef CENG_CARRAY2D_SAFE
		if ( _x < 0 || _y < 0 || _x >= myWidth || _y >= myHeight ) return myNullReference;
#else
		if ( _x < 0 ) _x = 0;
		if ( _y < 0 ) _y = 0;
		if ( _x >= myWidth ) _x = myWidth - 1;
		if ( _y >= myHeight ) _y = myHeight - 1;
#endif

		return myDataArray[ ( _y * myWidth ) + _x ];
	}


	inline const_reference At( int _x, int _y ) const
	{
#ifdef CENG_CARRAY2D_SAFE
		if ( _x < 0 || _y < 0 || _x >= myWidth || _y >= myHeight ) return myNullReference;
#else

		if ( _x < 0 ) _x = 0;
		if ( _y < 0 ) _y = 0;
		if ( _x >= myWidth ) _x = myWidth - 1;
		if ( _y >= myHeight ) _y = myHeight - 1;
#endif

		return myDataArray[ ( _y * myWidth ) + _x ];
	}



	inline reference Rand( int _x, int _y )
	{
#ifdef CENG_CARRAY2D_SAFE
		if ( _x < 0 || _y < 0 || _x >= myWidth || _y >= myHeight ) return myNullReference;
#endif
		return myDataArray[ ( _y * myWidth ) + _x ];
	}

	inline const_reference Rand( int _x, int _y ) const
	{
#ifdef CENG_CARRAY2D_SAFE
		if ( _x < 0 || _y < 0 || _x >= myWidth || _y >= myHeight ) return myNullReference;
#endif
		return myDataArray[ ( _y * myWidth ) + _x ];
	}

	void Rand( int _x, int _y, const _Ty& _who )
	{
		myDataArray[ ( _y * myWidth ) + _x ] = _who;
	}

	void Set( int _x, int _y, const _Ty& _who )
	{

		if ( _x > myWidth ) _x = myWidth;
		if ( _y > myHeight ) _y = myHeight;

		myDataArray[ ( _y * myWidth ) + _x ] = _who;
	}

	void Set( int _x, int _y, const CArray2D& _who )
	{
		int x, y;

		for ( y = 0; y <= myHeight; y++ )
		{
			for ( x = 0; x <= myWidth; x++ )
			{
				if ( x >= _x && x <= _x + _who.GetWidth() &&
					 y >= _y && y <= _y + _who.GetHeight() )
				{
					Set( x, y, _who.At( x - _x, y - _y ) );
				}
			}
		}
	}

	void Crop( const _Ty& _empty )
	{
		int left = myWidth;
		int right = 0;
		int top = myHeight;
		int bottom = 0;

		int x = 0;
		int y = 0;

		for ( y = 0; y <= myHeight; y++ )
		{
			for ( x = 0; x <= myWidth; x++ )
			{
				if ( At( x, y ) != _empty )
				{
					if ( x < left )		left	= x;
					if ( x > right )	right	= x;
					if ( y < top )		top		= y;
					if ( y > bottom )	bottom	= y;
				}
			}
		}

		Crop( left, top, right - left, bottom - top );
	}

	void Crop( int _x, int _y, int _w, int _h )
	{
		cassert(false);

	    /*
		std::vector< _Ty > tmpDataArray;

		tmpDataArray.resize( ( _w + 1 ) * ( _h + 1 ) );

		int x, y;
		for ( y = _y; y <= _y + _h; y++ )
		{
			for ( x = _x; x <= _x + _w; x++ )
			{
				tmpDataArray[ ( ( y - _y ) * _w ) + ( x - _x ) ] = At( x, y );
			}
		}

		myDataArray = tmpDataArray;
		myWidth = _w;
		myHeight = _h;*/
	}

	void Clear()
	{
		myWidth = 0;
		myHeight = 0;
		mySize = 0;
		myDataArray.clear();
	}

	bool Empty() const { return myDataArray.empty(); }

	CSafeArray< _Ty >& GetData() { return myDataArray; }
	const CSafeArray< _Ty >& GetData() const { return myDataArray; }

	CArray2D< _Ty, _A>* CopyCropped( int _x, int _y, int _w, int _h)
	{
		CArray2D< _Ty, _A>* result = new CArray2D< _Ty, _A >( _w, _h);

		int x, y;
		for ( y = _y; y < _y + _h; y++ )
		{
			for ( x = _x; x < _x + _w; x++ )
			{
				result->myDataArray[ ( ( y - _y ) * _w ) + ( x - _x ) ] = At( x, y );
			}
		}

		return result;
	}

private:

	void Allocate()
	{
		/*
		int n_size = (myWidth + 1) * ( myHeight + 1 );
		if( n_size != mySize )
		{
			mySize = n_size;
			myDataArray.resize( mySize + 1 );
		}*/

		int n_size = (myWidth) * ( myHeight );
		if( n_size != mySize )
		{
			mySize = n_size;
			myDataArray.resize( mySize + 1 );
		}
	}

	int myWidth;
	int myHeight;

	int mySize;

	CArray2DHelper	   myArraysLittleHelper;

	_Ty					myNullReference;

	CSafeArray< _Ty > myDataArray;
	// std::vector< _Ty > myDataArray;
};

//---------------------------------------------------------------------------------------------

#define DEF_RMask 0x000000FF
#define	DEF_GMask 0x0000FF00
#define	DEF_BMask 0x00FF0000
#define	DEF_AMask 0xFF000000
#define	DEF_RShift 0
#define	DEF_GShift 8
#define	DEF_BShift 16
#define	DEF_AShift 24


class CColorFloat
{
public:

	typedef unsigned int uint32;
	typedef unsigned char uint8;


	CColorFloat( float r = 0, float g = 0, float b = 0, float a = 1.f ) :
		r( r ),
		g( g ),
		b( b ),
		a( a ),
		multiplied_with_alpha( false )
	{ }

	CColorFloat( const CColorFloat& other ) : 
		r( other.r ),
		g( other.g ),
		b( other.b ),
		a( other.a ),
		multiplied_with_alpha( false )
	{ }


	CColorFloat( const uint32 clor ) :
		multiplied_with_alpha( false )
	{
		Set32( clor );	
	}
	
	//-------------------------------------------------------------------------
	// Operators
	CColorFloat operator+(const CColorFloat& other) const
	{
		return CColorFloat(r+other.r,g+other.g,b+other.b,a+other.a);
	}
	
	CColorFloat operator-(const CColorFloat& other) const
	{
		return CColorFloat(r-other.r,g-other.g,b-other.b,a-other.a);
	}
	
	CColorFloat operator*(const CColorFloat& other) const
	{
		return CColorFloat(r*other.r,g*other.g,b*other.b,a*other.a);
	}
	
	CColorFloat operator/(const CColorFloat& other) const
	{
		return CColorFloat(r/other.r,g/other.g,b/other.b,a/other.a);
	}
	
	CColorFloat operator*( float num ) const
	{
		return CColorFloat(r*num,g*num,b*num,a*num);
	}
	
	CColorFloat operator/( float  num ) const
	{
		return CColorFloat(r/num,g/num,b/num,a/num);
	}
	
	void operator+=(const CColorFloat& other)
	{
		r+=other.r;
		g+=other.g;
		b+=other.b;
		a+=other.a;
	}
	
	void operator-=(const CColorFloat& other)
	{
		r-=other.r;
		g-=other.g;
		b-=other.b;
		a-=other.a;
	}
	
	void operator*=( float  num)
	{
		r*=num;
		g*=num;
		b*=num;
		a*=num;
	}
	
	void operator/=(float num)
	{
		r/=num;
		g/=num;
		b/=num;
		a/=num;
	}
	
	bool FloatCompare( float v1, float v2, float e = 1.f / 512.f ) const
	{
		return ( fabs( v1 - v2 ) < e );

	}

	bool operator==(const CColorFloat& other) const
	{
		if (	FloatCompare( r, other.r ) &&
				FloatCompare( g, other.g ) &&
				FloatCompare( b, other.b ) &&
				FloatCompare( a, other.a ) ) return true;
		return false;
	}
	
	bool operator!=(const CColorFloat& other) const
	{
		return !operator==(other);
	}

	CColorFloat& operator=( const CColorFloat& other )
	{
		r = other.r;
		g = other.g;
		b = other.b;
		a = other.a;
		
		return *this;
	}

	float operator[]( int i ) const 
	{
		switch( i ) 
		{
		case 0:
			return GetR();
		case 1:
			return GetG();
		case 2:
			return GetB();
		case 3:
			return GetA();
		default:
			return 0;
		}
	}

	float& operator[]( int i )
	{
		static float dump_me = 0;
		switch( i ) 
		{
		case 0:
			return r;
		case 1:
			return g;
		case 2:
			return b;
		case 3:
			return a;
		default:
			return dump_me;
		}
	}

	///////////////////////////////////////////////////////////////////////////


	float GetR() const { return r; }
	float GetG() const { return g; }
	float GetB() const { return b; }
	float GetA() const { return a; }
	
	uint8 GetR8() const { return (uint8)(r * 255.0f); }
	uint8 GetG8() const { return (uint8)(g * 255.0f); }
	uint8 GetB8() const { return (uint8)(b * 255.0f); }
	uint8 GetA8() const { return (uint8)(a * 255.0f); }
	
	void SetR( float cr ) { r = cr; }
	void SetG( float cg ) { g = cg; }
	void SetB( float cb ) { b = cb; }
	void SetA( float ca ) { a = ca; }

	void Set8( uint8 r, uint8 g, uint8 b, uint8 a )
	{
		SetR( (float)(r / 255.f) );
		SetG( (float)(g / 255.f) );
		SetB( (float)(b / 255.f) );
		SetA( (float)(a / 255.f) );
	}

	void Set32( const uint32& color ) 
	{
		uint32 r32, g32, b32, a32;
		
		r32 = color & DEF_RMask;
		g32 = color & DEF_GMask;
		b32 = color & DEF_BMask;
		a32 = color & DEF_AMask;

		r = ( r32 >> DEF_RShift ) / 255.f;
		g = ( g32 >> DEF_GShift ) / 255.f;
		b = ( b32 >> DEF_BShift ) / 255.f;
		a = ( a32 >> DEF_AShift ) / 255.f;
	}

	// BUGBUG: Chechk should the return be r32 | g32 | b32 | a32 !?
	uint32 Get32() const
	{
		return ( GetR8() << DEF_RShift ) | 
			( GetG8() << DEF_GShift ) |
			( GetB8() << DEF_BShift ) |
			( GetA8() << DEF_AShift );
		/*
		uint32 r32, g32, b32, a32;
		
		r32 = GetR8() << DEF_RShift;
		g32 = GetG8() << DEF_GShift;
		b32 = GetB8() << DEF_BShift;
		a32 = GetA8() << DEF_AShift;

		return r32 | g32 | b32 | a32;*/
	}

	//=========================================================================

	CColorFloat GetMultipliedWithAlpha() const
	{
		CColorFloat result( *this );
		result.multiplied_with_alpha = true;
		
		/*if( this->multiplied_with_alpha )
			return result;*/

		result.r *= result.a;
		result.g *= result.a;
		result.b *= result.a;

		return result;
	}

	///////////////////////////////////////////////////////////////////////////
	float	r;	/*!<	The red component of this color		*/
	float	g;	/*!<	The green component of this color	*/
	float	b;	/*!<	The blue component of this color	*/
	float	a;	/*!<	The alpha component of this color	*/

	bool multiplied_with_alpha;
private:

	// static bool masks_initialized;

public:
	/*static uint32	DEF_RMask;
	static uint32	DEF_GMask;
	static uint32	DEF_BMask;
	static uint32	DEF_AMask;
	
	static uint8  DEF_RShift;
	static uint8  DEF_GShift;
	static uint8  DEF_BShift;
	static uint8  DEF_AShift;*/
	
};

template< typename T >
CColorFloat operator * ( T s, const CColorFloat& c) {
	return CColorFloat(s * c.r, s * c.g, s * c.b, s * c.a );
}

// -------------------------------- random ---------------------------------------
float Randomf( float low, float high );
int Random( int low, int high );

inline float Randomf( float low, float high )
{
	// REMEMBER TO INIT srand(time(NULL));
	return low+((high-low)*((float)std::rand() / (float)RAND_MAX) );
}


inline int Random( int low, int high )
{
	int t = high - low;
	return ( std::rand()%(t+1) ) + low;
}

} // end of namespace ceng

// ------ types -----
namespace types { 
	typedef ceng::math::CVector2< float >	vector2;
	typedef ceng::math::CVector2< double >	dvector2;
	typedef ceng::math::CVector2< int >		ivector2;

	typedef ceng::math::CVector2< int >		point;
	typedef ceng::CColorFloat				fcolor;

} // end of namespace types

#endif