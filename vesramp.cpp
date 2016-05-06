/*********
Building on vesicles.cpp 
This version makes temperature a local variable
to display a temperature ramp across the lattice.
************/

#include "includes.h"
#include <time.h>

// Influence list: t=0.0079sec/frame
int InfluenceList[8][14][2] =
{
	{ // ym = y-1
	  { -1, -2 },
	  { 0, -2 },
	  { 1, -2 },
	  {-1, -1},
	  { 0, -1},
	  { 1, -1},
	  { -1, 0 },
	  { 0, 0 },
	  { 1, 0 },
	  { -1, 1 },
	  { 0, 1 },
	  { 1, 1 }
	},
	{ // xm = x+1
	  { -1, -1 },
	  {  0, -1 },
	  {  1, -1 },
	  {  2, -1 },
	  { -1,  0 },
	  {  0,  0 },
	  {  1,  0 },
	  {  2,  0 },
	  { -1,  0 },
	  {  0,  1 },
	  {  1,  1 },
	  {  2,  1 }
	},
	{ // ym = y+1
	  { -1, -1 },
	  { 0, -1 },
	  { 1, -1 },
	  {-1, 0},
	  { 0, 0},
	  { 1, 0},
	  { -1, 1 },
	  { 0,  1 },
	  { 1,  1 },
	  { -1, 2 },
	  { 0,  2 },
	  { 1,  2 }
    },
	{ // xm = x-1
	  { -2, -1 },
	  { -1, -1 },
	  {  0, -1 },
	  {  1, -1 },
	  { -2,  0 },
	  { -1,  0 },
	  {  0,  0 },
	  {  1,  0 },
	  { -2,  0 },
	  { -1,  1 },
	  {  0,  1 },
	  {  1,  1 }
	},
	{ // xm = x-1, ym = y-1
	  { -2, -2 },
	  { -1, -2 },
	  {  0, -2 },
	  { -2, -1 },
	  { -1, -1 },
	  {  0, -1 },
	  {  1, -1  },
	  { -2, 0  },
	  { -1, 0  },
	  {  0, 0  },
	  {  1, 0  },
	  { -1, 1  },
	  { 0, 1  },
	  { 1, 1  }
	},
	{ // xm = x+1, ym = y-1
	  { 0, -2 },
	  { 1, -2 },
	  { 2, -2 },
	  { -1, -1 },
	  { 0, -1 },
	  { 1, -1 },
	  { 2, -1  },
	  { -1, 0  },
	  { 0, 0  },
	  { 1, 0  },
	  { 2, 0  },
	  { -1, 1  },
	  { 0, 1  },
	  { 1, 1  }
	},
	{ // xm = x-1, ym = y+1
	  { -1, -1 },
	  {  0, -1 },
	  {  1, -1 },
	  { -2, 0 },
	  { -1, 0 },
	  {  0, 0 },
	  {  1, 0  },
	  { -2, 1  },
	  { -1, 1  },
	  {  0, 1  },
	  {  1, 1  },
	  { -2, 2  },
	  { -1, 2  },
	  {  0, 2  }
    },
	{ // xm = x+1, ym = y+1
	  { -1, -1 },
	  { 0, -1 },
	  { 1, -1 },
	  { -1, 0 },
	  { 0, 0 },
	  { 1, 0 },
	  { 2, 0  },
	  { -1, 1  },
	  { 0, 1  },
	  { 1, 1  },
	  { 2, 1  },
	  { 0, 2  },
	  { 1, 2  },
	  { 2, 2  }
	},
};

class Cell
{
public:
    char phase;
    char bphase;
    float mx,my,bmx,bmy;
    float temp;
    void Backup();
    void Restore();
}; 

int XR,YR;
Cell *Grid;

double INTERFACE = 4.0;
double TENSION = 8.0;
double TEMP = 12.0;

double interaction[9] =
{
	-2,   -1,  1,
	-1,    0,  -1,
	1,    -1,  -2
};

unordered_map<string, double> pHash;

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

/* Uses exact patterns of neighborhoods 
 * t/frame = 0.0079sec
 * 
 * Does not work when energy involves non-linear terms in the neighbor counts and things like that ...
 */

double getHashedEnergy(int x, int y, int xm, int ym)
{
	char keyGrid[14];
	int xm2,ym2, cx,cy;
	int d;
	int ux,uy,vx,vy;
	string key = "";
	int i;
	
	if (abs(xm-x)+abs(ym-y)==1) // Straight
	{
		if (xm>x)
		{
			ux=1; uy=0;
			vx=0; vy=1;
		}
		else if (xm<x)
		{
			ux=-1; uy=0;
			vx=0; vy=-1;
		}
		else if (ym<y)
		{
			ux=0; uy=-1;
			vx=1; vy=0;
		}
		else
		{
			ux=0; uy=1;
			vx=-1; vy=0;
		}
		
		for (ym2=0;ym2<4;ym2++)
		{
			for (xm2=0;xm2<3;xm2++)
			{
				cx = MapX(x+(xm2-1)*vx + (ym2-2)*ux);
				cy = MapY(y+(xm2-1)*vy + (ym2-2)*uy);
								
				keyGrid[xm2+ym2*3] = Grid[cx+cy*XR].phase+'0';
			}
		}
		key = "S";
		for (i=0;i<12;i++)
			key = key + keyGrid[i];
	}
	else // Diagonal
	{
		if (xm>x)
		{
			if (ym>y)
			{
				cy=MapY(y);   cx=MapX(x+2); keyGrid[0]=Grid[cx+cy*XR].phase+'0';
				cy=MapY(y+1); cx=MapX(x+2); keyGrid[1]=Grid[cx+cy*XR].phase+'0';
				cy=MapY(y+2); cx=MapX(x+2); keyGrid[2]=Grid[cx+cy*XR].phase+'0';
				cy=MapY(y-1); cx=MapX(x+1); keyGrid[3]=Grid[cx+cy*XR].phase+'0';
				cy=MapY(y);   cx=MapX(x+1); keyGrid[4]=Grid[cx+cy*XR].phase+'0';
				cy=MapY(y+1); cx=MapX(x+1); keyGrid[5]=Grid[cx+cy*XR].phase+'0';
				cy=MapY(y+2); cx=MapX(x+1); keyGrid[6]=Grid[cx+cy*XR].phase+'0';
				cy=MapY(y-1); cx=MapX(x); keyGrid[7]=Grid[cx+cy*XR].phase+'0';
				cy=MapY(y);   cx=MapX(x); keyGrid[8]=Grid[cx+cy*XR].phase+'0';
				cy=MapY(y+1); cx=MapX(x); keyGrid[9]=Grid[cx+cy*XR].phase+'0';
				cy=MapY(y+2); cx=MapX(x); keyGrid[10]=Grid[cx+cy*XR].phase+'0';
				cy=MapY(y-1); cx=MapX(x-1); keyGrid[11]=Grid[cx+cy*XR].phase+'0';
				cy=MapY(y);   cx=MapX(x-1); keyGrid[12]=Grid[cx+cy*XR].phase+'0';
				cy=MapY(y+1); cx=MapX(x-1); keyGrid[13]=Grid[cx+cy*XR].phase+'0';				
			} else
			{
				cx=MapX(x);   cy=MapY(y-2); keyGrid[0]=Grid[cx+cy*XR].phase+'0';
				cx=MapX(x+1); cy=MapY(y-2); keyGrid[1]=Grid[cx+cy*XR].phase+'0';
				cx=MapX(x+2); cy=MapY(y-2); keyGrid[2]=Grid[cx+cy*XR].phase+'0';
				cx=MapX(x-1); cy=MapY(y-1); keyGrid[3]=Grid[cx+cy*XR].phase+'0';
				cx=MapX(x);   cy=MapY(y-1); keyGrid[4]=Grid[cx+cy*XR].phase+'0';
				cx=MapX(x+1); cy=MapY(y-1); keyGrid[5]=Grid[cx+cy*XR].phase+'0';
				cx=MapX(x+2); cy=MapY(y-1); keyGrid[6]=Grid[cx+cy*XR].phase+'0';
				cx=MapX(x-1); cy=MapY(y); keyGrid[7]=Grid[cx+cy*XR].phase+'0';
				cx=MapX(x);   cy=MapY(y); keyGrid[8]=Grid[cx+cy*XR].phase+'0';
				cx=MapX(x+1); cy=MapY(y); keyGrid[9]=Grid[cx+cy*XR].phase+'0';
				cx=MapX(x+2); cy=MapY(y); keyGrid[10]=Grid[cx+cy*XR].phase+'0';
				cx=MapX(x-1); cy=MapY(y+1); keyGrid[11]=Grid[cx+cy*XR].phase+'0';
				cx=MapX(x);   cy=MapY(y+1); keyGrid[12]=Grid[cx+cy*XR].phase+'0';
				cx=MapX(x+1); cy=MapY(y+1); keyGrid[13]=Grid[cx+cy*XR].phase+'0';
			}
		}
		else 
		{
			if (ym>y)
			{
				cx=MapX(x);   cy=MapY(y+2); keyGrid[0]=Grid[cx+cy*XR].phase+'0';
				cx=MapX(x-1); cy=MapY(y+2); keyGrid[1]=Grid[cx+cy*XR].phase+'0';
				cx=MapX(x-2); cy=MapY(y+2); keyGrid[2]=Grid[cx+cy*XR].phase+'0';
				cx=MapX(x+1); cy=MapY(y+1); keyGrid[3]=Grid[cx+cy*XR].phase+'0';
				cx=MapX(x);   cy=MapY(y+1); keyGrid[4]=Grid[cx+cy*XR].phase+'0';
				cx=MapX(x-1); cy=MapY(y+1); keyGrid[5]=Grid[cx+cy*XR].phase+'0';
				cx=MapX(x-2); cy=MapY(y+1); keyGrid[6]=Grid[cx+cy*XR].phase+'0';
				cx=MapX(x+1); cy=MapY(y);   keyGrid[7]=Grid[cx+cy*XR].phase+'0';
				cx=MapX(x);   cy=MapY(y);   keyGrid[8]=Grid[cx+cy*XR].phase+'0';
				cx=MapX(x-1); cy=MapY(y);   keyGrid[9]=Grid[cx+cy*XR].phase+'0';
				cx=MapX(x-2); cy=MapY(y);   keyGrid[10]=Grid[cx+cy*XR].phase+'0';
				cx=MapX(x+1); cy=MapY(y-1); keyGrid[11]=Grid[cx+cy*XR].phase+'0';
				cx=MapX(x);   cy=MapY(y-1); keyGrid[12]=Grid[cx+cy*XR].phase+'0';
				cx=MapX(x-1); cy=MapY(y-1); keyGrid[13]=Grid[cx+cy*XR].phase+'0';
			} else
			{
				cy=MapY(y);   cx=MapX(x-2); keyGrid[0]=Grid[cx+cy*XR].phase+'0';
				cy=MapY(y-1); cx=MapX(x-2); keyGrid[1]=Grid[cx+cy*XR].phase+'0';
				cy=MapY(y-2); cx=MapX(x-2); keyGrid[2]=Grid[cx+cy*XR].phase+'0';
				cy=MapY(y+1); cx=MapX(x-1); keyGrid[3]=Grid[cx+cy*XR].phase+'0';
				cy=MapY(y);   cx=MapX(x-1); keyGrid[4]=Grid[cx+cy*XR].phase+'0';
				cy=MapY(y-1); cx=MapX(x-1); keyGrid[5]=Grid[cx+cy*XR].phase+'0';
				cy=MapY(y-2); cx=MapX(x-1); keyGrid[6]=Grid[cx+cy*XR].phase+'0';
				cy=MapY(y+1); cx=MapX(x);   keyGrid[7]=Grid[cx+cy*XR].phase+'0';
				cy=MapY(y);   cx=MapX(x);   keyGrid[8]=Grid[cx+cy*XR].phase+'0';
				cy=MapY(y-1); cx=MapX(x);   keyGrid[9]=Grid[cx+cy*XR].phase+'0';
				cy=MapY(y-2); cx=MapX(x);   keyGrid[10]=Grid[cx+cy*XR].phase+'0';
				cy=MapY(y+1); cx=MapX(x+1); keyGrid[11]=Grid[cx+cy*XR].phase+'0';
				cy=MapY(y);   cx=MapX(x+1); keyGrid[12]=Grid[cx+cy*XR].phase+'0';
				cy=MapY(y-1); cx=MapX(x+1); keyGrid[13]=Grid[cx+cy*XR].phase+'0';				
			}
		}
		
		key = "D";
		for (i=0;i<14;i++)
			key = key + keyGrid[i];
	}
	
	char other = Grid[xm+ym*XR].phase, local = Grid[x+y*XR].phase;
	Grid[x+y*XR].Backup(); Grid[xm+ym*XR].Backup();
	
	if (pHash.count(key))
	{
		Grid[x+y*XR].phase = other; Grid[xm+ym*XR].phase = local;
		return pHash[key];
	}
	else
	{
		float E0=0,E1=0;
		int xm2,ym2;
		
		for (ym2=y-2;ym2<=y+2;ym2++)
			for (xm2=x-2;xm2<=x+2;xm2++)
			{
				E0+=localEnergy(MapX(xm2),MapY(ym2));
			}

		Grid[x+y*XR].phase = other; Grid[xm+ym*XR].phase = local;

		for (ym2=y-2;ym2<=y+2;ym2++)
			for (xm2=x-2;xm2<=x+2;xm2++)
			{
				E1+=localEnergy(MapX(xm2),MapY(ym2));
			}
			
		pHash[key]=E1-E0;
		return E1-E0;
	}
}

/* Uses sum of neighbors */
// This doesn't work with all Hamiltonians that only use local information!!!
/* t = 0.00128 sec/frame */
/*double getHashedEnergy(int x, int y, int xm, int ym)
{
	char keyGrid[14];
	int xm2,ym2, cx,cy;
	int d;
	int ux,uy,vx,vy;
	char noil=0, nlipid=0;
	int key;
	int i;
	int mult = 1;
	
	if (abs(xm-x)+abs(ym-y)==1) key = key + 0;
	else key = key + 1;
	mult*=2;
	
	key += Grid[x+y*XR].phase*mult; mult*=3;
	key += Grid[xm+ym*XR].phase*mult; mult*=3;
	
	for (ym2=y-1;ym2<=y+1;ym2++)
		for (xm2=x-1;xm2<=x+1;xm2++)
		{
			cx = MapX(xm2); cy=MapY(ym2);
			char val = Grid[cx+cy*XR].phase;
			
			if ((xm2!=x)||(ym2!=y))
			{
				if (val==0) noil++;
				else if (val==1) nlipid++;
			}
		}
	key = key + noil*mult; mult*=9;
	key = key + nlipid*mult; mult*=9;
	
	noil = 0; nlipid = 0;
	for (ym2=ym-1;ym2<=ym+1;ym2++)
		for (xm2=xm-1;xm2<=xm+1;xm2++)
		{
			cx = MapX(xm2); cy=MapY(ym2);
			char val = Grid[cx+cy*XR].phase;
			
			if ((xm2!=x)||(ym2!=y))
			{
				if (val==0) noil++;
				else if (val==1) nlipid++;
			}
		}
	key = key + noil*mult; mult*=9;
	key = key + nlipid*mult; mult*=9;
	
	char other = Grid[xm+ym*XR].phase, local = Grid[x+y*XR].phase;
	Grid[x+y*XR].Backup(); Grid[xm+ym*XR].Backup();
	
	if (pHash.count(key))
	{
		Grid[x+y*XR].phase = other; Grid[xm+ym*XR].phase = local;
		return pHash[key];
	}
	else
	{
		float E0=0,E1=0;
		int xm2,ym2;
		
		for (ym2=y-2;ym2<=y+2;ym2++)
			for (xm2=x-2;xm2<=x+2;xm2++)
			{
				E0+=localEnergy(MapX(xm2),MapY(ym2));
			}

		Grid[x+y*XR].phase = other; Grid[xm+ym*XR].phase = local;

		for (ym2=y-2;ym2<=y+2;ym2++)
			for (xm2=x-2;xm2<=x+2;xm2++)
			{
				E1+=localEnergy(MapX(xm2),MapY(ym2));
			}
			
		pHash[key]=E1-E0;
		return E1-E0;
	}
}*/

void Iterate()
{
	int x,y;
	int xm,ym,i;
	
	x=rand()%XR; y=rand()%YR;
	
	int stop=0;
	int d = rand()%8;
	
	switch (d)
	{
		case 0: ym=y-1; xm=x; break;
		case 1: ym=y; xm=x+1; break;
		case 2: ym=y+1; xm=x; break;
		case 3: ym=y; xm=x-1; break;
		case 4: xm=x-1; ym=y-1; break;
		case 5: xm=x+1; ym=y-1; break;
		case 6: xm=x-1; ym=y+1; break;
		case 7: xm=x+1; ym=y+1; break;
	}

	xm=MapX(xm); ym=MapY(ym);
	
	int local, other;
	
	local=Grid[x+y*XR].phase; other=Grid[xm+ym*XR].phase;
	
	if (local!=other)
	{
		/*
		double dE;
		int xm2,ym2;
		
		dE = getHashedEnergy(x,y,xm,ym);

		if (dE>0)
		{
			if (rand()%1000001 >= 1000000.0 * exp(-dE/TEMP))
			{
				Grid[x+y*XR].Restore();
				Grid[xm+ym*XR].Restore();
			}
		}
		*/	
		
		float E0=0,E1=0;
		int xm2,ym2;
		int nInf;
		
		if (d<4) nInf = 12; else nInf = 14;
		
		for (i=0;i<nInf;i++)
		{
			E0+=localEnergy(MapX(x+InfluenceList[d][i][0]),MapY(y+InfluenceList[d][i][1]));
		}

		Grid[x+y*XR].Backup(); Grid[xm+ym*XR].Backup();
		Grid[x+y*XR].phase = other; Grid[xm+ym*XR].phase = local;

		for (i=0;i<nInf;i++)
		{
			E1+=localEnergy(MapX(x+InfluenceList[d][i][0]),MapY(y+InfluenceList[d][i][1]));
		}
			
		float dE=E1-E0;
			
		if (dE>0)
		{
			if (rand()%1000001 >= 1000000.0 * exp(-dE/Grid[x+y*XR].temp))
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
			
//			ScreenBuf[(xm+ym*XRes)*Bpp]=b;
//			ScreenBuf[(xm+ym*XRes)*Bpp+1]=g;
//			ScreenBuf[(xm+ym*XRes)*Bpp+2]=r;
			ScreenBuf[(xm+ym*XRes)*Bpp+3]=b;
			ScreenBuf[(xm+ym*XRes)*Bpp+2]=g;
			ScreenBuf[(xm+ym*XRes)*Bpp+1]=r;
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
			Grid[x+y*XR].temp=(TEMP/2) + x*TEMP/(float)(XR);
        }            
	for (y=0;y<YR;y++)
		for (x=0;x<XR;x++)
		{
			int r=sqrt(pow(x-XR/2,2)+pow(y-YR/2,2));
			
//			if (r<30) Grid[x+y*XR].phase = 0;
//			else Grid[x+y*XR].phase=2;
			Grid[x+y*XR].phase=2;

			if (rand()%5==0) Grid[x+y*XR].phase=0;
			if (rand()%10==0) Grid[x+y*XR].phase=1;
		}
}

int main(int argc, char **argv)
{
	int t1,t2;
	
	Init();
	XRes=XR*4; YRes=YR*4;
	Bpp=4;
	ScreenBuf=(unsigned char*)malloc(XRes*YRes*Bpp);
	InitSDL();
	
	while (1)
	{
		int Ch=ReadKey();
		
		if (Ch=='q') return 0;
		
//		t1=clock();
		
		for (int i=0;i<XR*YR*10;i++)
			Iterate();
			
//		t2=clock();
		
		//printf("%.6g\n",(t2-t1)/(100.0*CLOCKS_PER_SEC));
		
//		TEMP*=0.999; if (TEMP<8) TEMP=8;
		Render();
	}
}
