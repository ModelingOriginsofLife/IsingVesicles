#include "includes.h"

class Cell
{
	public:
		char phase;
		char bphase;
		float mx,my,bmx,bmy;
		
		void Backup();
		void Restore();
}; 

int XR,YR;
Cell *Grid;

double INTERFACE = 4.0;
double TENSION = 8.0;
double TEMP = 8.0;

double interaction[9] =
{
	-2,   -1,  1,
	-1,    0,  -1,
	1,    -1,  -2
};

void Cell::Backup()
{
	bphase = phase; bmx=mx; bmy=my;
}

void Cell::Restore()
{
	phase = bphase; mx=bmx; my=bmy;
}

int MapX(int xm)
{
	if (xm<0) xm+=XR; if (xm>=XR) xm-=XR;
	
	return xm;
}

int MapY(int ym)
{
	if (ym<0) ym+=YR; if (ym>=YR) ym-=YR;
	
	return ym;
}

float localEnergy(int x, int y)
{
	int xm,ym;
	float E=0;
	int local = Grid[x+y*XR].phase;
	int neighbor, bdylayer = 0, noil = 0, nwater = 0;
	
	for (ym=y-1;ym<=y+1;ym++)
		for (xm=x-1;xm<=x+1;xm++)
		{			
			if ((xm!=x)||(ym!=y))
			{
				int xm2=xm,ym2=ym;
				
				xm2=MapX(xm2); ym2=MapY(ym2);
				
				neighbor = Grid[xm2+ym2*XR].phase;
				
				E+=interaction[local+neighbor*3];
				
				if (neighbor == 0) noil++;
				if (neighbor == 1) bdylayer++;
				if (neighbor == 2) nwater++;
			}
		}
		
	if (local == 1)
	{
		E += TENSION * (bdylayer - 2)*(bdylayer - 2) + INTERFACE * fabs(noil-nwater);
	}
		
	return E;
}

void Iterate()
{
	int x,y;
	int xm,ym;
	
	x=rand()%XR; y=rand()%YR;
	
	int stop=0;
	float dp;
	int xm0,ym0;
	
	xm=(x+rand()%3-1); ym=(y+rand()%3-1);
	xm0=xm; ym0=ym;
	xm=MapX(xm); ym=MapY(ym);
	
	int local, other;
	
	local=Grid[x+y*XR].phase; other=Grid[xm+ym*XR].phase;
	
	if (local!=other)
	{
		float E0=0, E1=0, othermx, othermy;
		int xm2,ym2;
		
		for (ym2=y-2;ym2<=y+2;ym2++)
			for (xm2=x-2;xm2<=x+2;xm2++)
			{
				E0+=localEnergy(MapX(xm2),MapY(ym2));
			}
			
		Grid[x+y*XR].Backup(); Grid[xm+ym*XR].Backup();
		
		Grid[x+y*XR].phase = other; Grid[xm+ym*XR].phase = local;
		
		for (ym2=y-2;ym2<=y+2;ym2++)
			for (xm2=x-2;xm2<=x+2;xm2++)
			{
				E1+=localEnergy(MapX(xm2),MapY(ym2));
			}
		
		if (E1>E0)
		{
			if (rand()%1000001 >= 1000000 * exp(-(E1-E0)/TEMP))
			{
				Grid[x+y*XR].Restore();
				Grid[xm+ym*XR].Restore();
			}		
		}
	}
}

void Render()
{
	int x,y,xm,ym;
	int r,g,b;
	
	for (y=0;y<YR;y++)
		for (x=0;x<XR;x++)
		{
			for (ym=4*y;ym<4*(y+1);ym++)
				for (xm=4*x;xm<4*(x+1);xm++)
				{
			switch (Grid[x+y*XR].phase)
			{
				case 0: r=255; g=0; b=0; break;
				case 1: r=0; g=255; b=0; break;
				case 2: r=0; g=0; b=255; break;
			}
			
			ScreenBuf[(xm+ym*XRes)*Bpp]=b;
			ScreenBuf[(xm+ym*XRes)*Bpp+1]=g;
			ScreenBuf[(xm+ym*XRes)*Bpp+2]=r;
		}
		}
		
	BlitBuf(ScreenBuf,0,0,XRes,YRes);
}

int getRandomPhase()
{
	int r = rand()%100;
	
	if (r<20) return 1;
	else if (r<45) return 0;
	else return 2;
}

void Init()
{
	XR = 128; YR = 128;
	
	Grid=(Cell*)malloc(XR*YR*sizeof(Cell));
	
	int x,y;
	
	for (y=0;y<YR;y++)
		for (x=0;x<XR;x++)
		{
			int r=sqrt(pow(x-XR/2,2)+pow(y-YR/2,2));
			
			if (r<30) Grid[x+y*XR].phase = 0;
			else Grid[x+y*XR].phase=2;
			
			if (rand()%20==0) Grid[x+y*XR].phase=1;
		}
}

int main(int argc, char **argv)
{
	Init();
	XRes=XR*4; YRes=YR*4;
	Bpp=4;
	ScreenBuf=(unsigned char*)malloc(XRes*YRes*Bpp);
	InitSDL();
	
	while (1)
	{
		int Ch=ReadKey();
		
		if (Ch=='q') return 0;
		
		for (int i=0;i<XR*YR*10;i++)
			Iterate();
			
		Render();
	}
}
