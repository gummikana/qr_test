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


void QRTest::Init()
{
	DefaultApplication::Init();
	Poro()->GetGraphics()->SetFillColor( poro::GetFColor( 248.f / 255.f, 245.f / 255.f, 236.f / 255.f, 1.f ) );

	mSpriteContainer = new as::Sprite;
	mDebugLayer.reset( new DebugLayer );

	// t4, t5, t12, t7
	DoCard( "test/t7.jpg", "cards/h12.png", "test_image.jpg" );

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
	if( config_display_wireframe && mDebugMarkers.empty() == false )
	{
		// int work_width = 1481;

		float scale = ( 1024.f / 1481.f );

		types::vector2 center_p(0,0);
		for( std::size_t k = 0; k < mDebugMarkers.size(); ++k )	
		{
			center_p += mDebugMarkers[k].center;
			int j, i;
			for (j=4-1,i=0; i < 4; j=i++) 
			{
				DrawLine( graphics, 
					scale * ( mDebugMarkers[k].corners[i] ),  
					scale * ( mDebugMarkers[k].corners[j] ),
					poro::GetFColor( 0, 1, 0, 0.5f ),
					NULL );
			}
		}


		/*
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

			*/
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

