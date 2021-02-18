#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <limits>


const double pi = 3.14159265359;

void normalize(double n[3])
{
	double l;
	l = std::sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);
	if( l )
	{
		n[0] /= l;
		n[1] /= l;
		n[2] /= l;
	}
}

//~ const double noind = std::numeric_limits<double>::quiet_NaN();
const double noind = 0.0;
int N = 128;
#define ind(r,c) ((r)*N + (c))


void rand2d(double * arr, int N, int Nl, int Nr, int Nb, int Nt, double t)
{
	//std::cout << "Nl:Nr " << Nl <<":" << Nr << " Nb:Nt " << Nb << ":" << Nt << std::endl;
	if(  Nr - Nl < 2 && Nt - Nb < 2 ) 
	{
		//std::cout << "exit" << std::endl;
		return;
	}
	//const double t = 0.15;
	int Nk = (Nb+Nt)/2;
	int Nm = (Nl+Nr)/2;
	//std::cout << "Nk " << Nk << " Nm " << Nm << std::endl;
	double lb = arr[ind(Nb,Nl)];
	double rb = arr[ind(Nb,Nr)];
	double lt = arr[ind(Nt,Nl)];
	double rt = arr[ind(Nt,Nr)];
	if( lb != lb || rb != rb || lt != lt || rt != rt ) throw -1;
	if( arr[ind(Nk,Nl)] != arr[ind(Nk,Nl)] ) arr[ind(Nk,Nl)] = 0.5*(lb + lt) + (2*(rand()*1.0/RAND_MAX)-1)*t;
	if( arr[ind(Nk,Nr)] != arr[ind(Nk,Nr)] ) arr[ind(Nk,Nr)] = 0.5*(rb + rt) + (2*(rand()*1.0/RAND_MAX)-1)*t;
	if( arr[ind(Nb,Nm)] != arr[ind(Nb,Nm)] ) arr[ind(Nb,Nm)] = 0.5*(lb + rb) + (2*(rand()*1.0/RAND_MAX)-1)*t;
	if( arr[ind(Nt,Nm)] != arr[ind(Nt,Nm)] ) arr[ind(Nt,Nm)] = 0.5*(lt + rt) + (2*(rand()*1.0/RAND_MAX)-1)*t;
	arr[ind(Nk,Nm)] = 0.25*(lb+rb+lt+rt) + (2*(rand()*1.0/RAND_MAX)-1)*t;
	rand2d(arr,N,Nl,Nm,Nb,Nk,t*0.5); rand2d(arr,N,Nm,Nr,Nb,Nk,t*0.5);
	rand2d(arr,N,Nl,Nm,Nk,Nt,t*0.5); rand2d(arr,N,Nm,Nr,Nk,Nt,t*0.5);
}

void init2d(double * arr, int N, double mint, double maxt)
{
	for(int k = 0; k < N*N; ++k) arr[k] = noind;
	arr[ind(0  ,0  )] = mint+(rand()*1.0/RAND_MAX)*(maxt-mint);
	arr[ind(0  ,N-1)] = mint+(rand()*1.0/RAND_MAX)*(maxt-mint);
	arr[ind(N-1,0  )] = mint+(rand()*1.0/RAND_MAX)*(maxt-mint);
	arr[ind(N-1,N-1)] = mint+(rand()*1.0/RAND_MAX)*(maxt-mint);
}

double intrp2d(const double * arr, int N, double x, double y)
{
	int n = ceil(x*(N-1));
	int m = ceil(y*(N-1));
	if( n == 0 ) n = 1;
	if( m == 0 ) m = 1;
	double dh = 1.0/(double)(N-1);
	double kx = (x-(n-1)*dh)/dh;
	double ky = (y-(m-1)*dh)/dh;
	//if( kx < 0 || kx > 1 ) std::cout << "bad kx: " << kx << " x is " << x << " n is " << n << " dh is " << dh << " N is " << N << std::endl;
	//if( ky < 0 || ky > 1 ) std::cout << "bad ky: " << ky << " y is " << y << " m is " << m << " dh is " << dh << " N is " << N << std::endl;
	double lb = arr[ind(m-1,n-1)];
	double rb = arr[ind(m-1,n+0)];
	double lt = arr[ind(m+0,n-1)];
	double rt = arr[ind(m+0,n+0)];
	if( lb != lb || rb != rb || lt != lt || rt != rt ) throw -1;
	return (1-ky)*(lb*(1-kx) + rb*kx) + ky*(lt*(1-kx) + rt*kx);
}

double intrp2dx(const double * arr, int N, double x, double y)
{
	int n = ceil(x*(N-1));
	int m = ceil(y*(N-1));
	if( n == 0 ) n = 1;
	if( m == 0 ) m = 1;
	double dh = 1.0/(double)(N-1);
	double lb = arr[ind(m-1,n-1)];
	double rb = arr[ind(m-1,n+0)];
	double lt = arr[ind(m+0,n-1)];
	double rt = arr[ind(m+0,n+0)];
	return 0.5*((rt+rb) - (lt+lb))/dh;
}

double intrp2dy(const double * arr, int N, double x, double y)
{
	int n = ceil(x*(N-1));
	int m = ceil(y*(N-1));
	if( n == 0 ) n = 1;
	if( m == 0 ) m = 1;
	double dh = 1.0/(double)(N-1);
	double lb = arr[ind(m-1,n-1)];
	double rb = arr[ind(m-1,n+0)];
	double lt = arr[ind(m+0,n-1)];
	double rt = arr[ind(m+0,n+0)];
	return 0.5*((rt+lt) - (rb+lb))/dh;
}


void transform(const std::vector<double> & map, double xyz[3], const double max[3], const double min[3], double & ztop, double & zbottom, double nrmtop[3], double nrmbottom[3])
{
	double x = (xyz[0]-min[0])/(max[0]-min[0]), y = (xyz[1]-min[1])/(max[1]-min[1]), z = (xyz[2]-min[2])/(max[2]-min[2]);
	if( x < 0 || x > 1 ) {throw -1; std::cout << "x: " << x << " xyz " << xyz[0] << " " << xyz[1] << " " << xyz[2] << " unit " << x << " " << y << " " << z << " min " << min[0] << " " << min[1] << " " << min[2] << " max " << max[0] << " " << max[1] << " " << max[2] << std::endl;}
	if( y < 0 || y > 1 ) std::cout << "y: " << y << " xyz " << xyz[0] << " " << xyz[1] << " " << xyz[2] << " unit " << x << " " << y << " " << z << " min " << min[0] << " " << min[1] << " " << min[2] << " max " << max[0] << " " << max[1] << " " << max[2] << std::endl;
	if( z < 0 || z > 1 ) std::cout << "z: " << z << " xyz " << xyz[0] << " " << xyz[1] << " " << xyz[2] << " unit " << x << " " << y << " " << z << " min " << min[0] << " " << min[1] << " " << min[2] << " max " << max[0] << " " << max[1] << " " << max[2] << std::endl;
	double dztopdx = 0, dztopdy = 0;
	double dzbottomdx = 0, dzbottomdy = 0;
	double shift = 0;
	//zbottom = x*x*4-y*y*4-sin(6*x)*8 - cos(4*y)*4 - x*15;
	//ztop = x*x*4-y*y*4-sin(6*x)*8 - cos(4*y)*4 - y*15 + 15;
	
	
	//zbottom = std::sin(y*pi *2)*0.2;// ((x - 0.5)*(x - 0.5) + (y - 0.5)*(y - 0.5))*0.4;
	//ztop = std::sin(y*pi *2)*0.2+1;// 1 + cos(4 * y)*0.2 + ((x - 0.5)*(x - 0.5) + (y - 0.5)*(y - 0.5))*0.4 + 1 + cos(4 * y)*0.5;
	//dzbottomdx = 0;
	//dzbottomdy = 2*pi*std::cos(y*pi*2)*0.2;
	//dztopdx = 0;
	//dztopdy = 2*pi*std::cos(y*pi*2)*0.2;
	
	shift = intrp2d(&map[0],N,x,y);
	zbottom = shift;
	ztop = 1+shift;
	dzbottomdx = dztopdx = intrp2dx(&map[0],N,x,y);
	dzbottomdy = dztopdy = intrp2dy(&map[0],N,x,y);
	
	
	
	//zbottom = 0;// ((x-0.5)*(x-0.5) + (y-0.5)*(y-0.5))*0.4;
	//ztop = 1 + 1*y;//cos(4*y)*0.2;// + ((x-0.5)*(x-0.5) + (y-0.5)*(y-0.5))*0.4 + 1 + cos(4*y)*0.5;
	
	//zbottom = 0;// ((x-0.5)*(x-0.5) + (y-0.5)*(y-0.5))*0.4;
	//ztop = 1 + cos(4*y)*0.2;// + ((x-0.5)*(x-0.5) + (y-0.5)*(y-0.5))*0.4 + 1 + cos(4*y)*0.5;
	
	
	
	
	
	nrmtop[0] = dztopdx;
	nrmtop[1] = dztopdy;
	nrmtop[2] = -1;
	normalize(nrmtop);
	nrmbottom[0] = dzbottomdx;
	nrmbottom[1] = dzbottomdy;
	nrmbottom[2] = -1;
	normalize(nrmbottom);
}

void changexyz(double xyz[3], double max[3], double min[3], double ztop, double zbottom, double xyzout[3])
{
	double tz = (xyz[2]-min[2])/(max[2]-min[2]), tznew;
	tznew = (ztop-zbottom)*tz + zbottom;
	xyzout[0] = xyz[0];
	xyzout[1] = xyz[1];
	xyzout[2] = tznew*(max[2]-min[2]) + min[2];
}


struct quat
{
	double vec[3];
	double w;
};

double quatnorm(const quat & a)
{
	return a.vec[0]*a.vec[0] + a.vec[1]*a.vec[1]
		 + a.vec[2]*a.vec[2] + a.w*a.w;
}

quat quatconj(const quat & a)
{
	quat ret = a;
	ret.vec[0] = -ret.vec[0];
	ret.vec[1] = -ret.vec[1];
	ret.vec[2] = -ret.vec[2];
	return ret;
}

void transposemxn(const double * a, double * b, int m, int n)
{
	for(int j = 0; j < m; j++)
	{
		for(int i = 0; i < n; i++)
			b[i*m+j] = a[j*n+i];
	}
}

void matmulmxnxk(const double * a, const double * b, double * out, int m, int n, int k)
{
	int i,j,l;
	for(i = 0; i < n; i++)
		for(j = 0; j < m; j++)
		{
			out[j*n+i] = 0.0;
			for(l = 0; l < k; l++)
				out[j*n+i] += a[j*k+l]*b[l*n+i];
		}
}

void rotate_tensor(double nrm[3], const double Kin[3], double Kout[6])
{
	double qmat[9], qmat_t[9], prod[9], qnrm;
	struct quat q, qr;
	
	q.vec[0] =-nrm[1];
	q.vec[1] = nrm[0];
	q.vec[2] = 0;
	q.w = 1.0 - nrm[2];

	qnrm = quatnorm(q);
	q.vec[0] /= qnrm;
	q.vec[1] /= qnrm;
	q.w /= qnrm;
	qr = quatconj(q);

	qmat[0] = (q.w*q.w + q.vec[0]*q.vec[0] - q.vec[1]*q.vec[1]) *qnrm;
	qmat[1] = 2.*(q.vec[0]*q.vec[1] )                           *qnrm;
	qmat[2] = 2.*(q.w*q.vec[1])                                 *qnrm;
	
	qmat[3] = 2.*(q.vec[0]*q.vec[1])                            *qnrm;
	qmat[4] = (q.w*q.w - q.vec[0]*q.vec[0] + q.vec[1]*q.vec[1]) *qnrm;
	qmat[5] = 2.*( - q.w*q.vec[0])                              *qnrm;
	
	qmat[6] = 2.*( - q.w*q.vec[1])                              *qnrm;
	qmat[7] = 2.*( + q.w*q.vec[0])                              *qnrm;
	qmat[8] = (q.w*q.w - q.vec[0]*q.vec[0] - q.vec[1]*q.vec[1]) *qnrm;
	
	transposemxn(qmat,qmat_t,3,3);
	
	qmat[0] *= Kin[0];
	qmat[1] *= Kin[0];
	qmat[2] *= Kin[0];
	
	qmat[3] *= Kin[1];
	qmat[4] *= Kin[1];
	qmat[5] *= Kin[1];
	
	qmat[6] *= Kin[2];
	qmat[7] *= Kin[2];
	qmat[8] *= Kin[2];

	matmulmxnxk(qmat_t,qmat,prod,3,3,3);

	Kout[0] = prod[0];
	Kout[1] = prod[1];
	Kout[2] = prod[2];
	Kout[3] = prod[4];
	Kout[4] = prod[5];
	Kout[5] = prod[8];
}
/*
void write_points_vtk(int lnx, int rnx, int refx,
					  int lny, int rny, int refy,
					  int lnz, int rnz, int refz,
					  double max[3], double min[3])
{
	double ztop, double zbottom, xyzout[3], xyz[3], nrmtop[3], nrmbottom[3];
	for(int k = lnz; k <= rnz; ++k)
	{
		for(int kr = 0; kr < (k < rnz ? refz : 1); ++kr)
		{
			for(int j = lny; j <= rny; ++j)
			{
				for(int jr = 0; jr < (j < rny ? refy : 1); jr++)
				{
					for(int i = lnx; i <= rnx; ++i)
					{
						for(int ir = 0; ir < (i < rnx ? refx : 1); ir++)
						{
							//bottom point
							xyz[0] = 240.0 * (i * refx + ir) / ( 60.0 * refx);
							xyz[1] = 440.0 * (j * refy + jr) / (220.0 * refy);
							xyz[2] = 340.0 * (k * refz + kr) / ( 85.0 * refz);
							transform(xyz,max,min,ztop,zbottom,nrmtop,nrmbottom);
							changexyz(xyz,max,min,ztop,zbottom,xyzout);
							fvtk << xyzout[0] << " " << xyzout[1] << " " << xyzout[2] << std::endl;
						}
					}
				}
			}
		}
	}
}
*/
int main(int argc, char *argv[]) 
{
	if (argc < 2)
	{
		std::cout << "Usage: " << argv[0] ;
		std::cout << " output.grdecl [deformation=0.5] [lnx=0 rnx=60 lny=0 rny=220 lnz=0 rnz=85] [refine_x=1] [refine_y=1] [refine_z=1] [write_vtk=1(1:hex,2:tet,3:hex(vtu),4:tet(vtu))]" << std::endl;
		return -1;
	}

	int wvtk = 1;
	std::cout << "Opening " << argv[1] << " for output." << std::endl;
	std::ofstream f(argv[1]), fvtk;
	
	if( f.fail() )
	{
		std::cout << "Cannot open " << argv[1] << " for writing!" << std::endl;
		return -1;
	}
	
	
	
	double value, ztop, zbottom, xyz[3], xyzout[3], nrmtop[3], nrmbottom[3], nrm[3], Kin[3], Kout[6];
	double max[3] = {240,440,340}, min[3] = {0,0,0};
	double deformation = 0.5;
	int nx = 60, ny = 220, nz = 85, m;
	int lnx = 0, lny = 0, lnz = 0;
	int rnx = nx, rny = ny, rnz = nz;
	int refx = 1, refy = 1, refz = 1;
	size_t nout = 0;
	bool vtu = false;
	
	if( argc > 2 ) deformation = atof(argv[2]);
	
	std::vector<double> map(N*N,0.0);
	init2d(&map[0],N,0.0,deformation);
	rand2d(&map[0],N,0,N-1,0,N-1,deformation*0.5);
	
	if( argc > 3  ) lnx  = atoi(argv[3]);
	if( argc > 4  ) rnx  = atoi(argv[4]);
	if( argc > 5  ) lny  = atoi(argv[5]);
	if( argc > 6  ) rny  = atoi(argv[6]);
	if( argc > 7  ) lnz  = atoi(argv[7]);
	if( argc > 8  ) rnz  = atoi(argv[8]);
	if( argc > 9  ) refx = atoi(argv[9]);
	if( argc > 10 ) refy = atoi(argv[10]);
	if( argc > 11 ) refz = atoi(argv[11]);
	if( argc > 12 ) wvtk = atoi(argv[12]);
	
	if( wvtk )
	{
		std::string ext = "vtk";
		if( wvtk == 3 || wvtk == 4 ) 
		{
			ext = "vtu";
			vtu = true;
		}
		std::cout << "Opening " << argv[1] << "." << ext << " for output." << std::endl;
		if( wvtk == 1 || wvtk == 3)
			std::cout << "Writing hexahedral mesh." << std::endl;
		else if( wvtk == 2 || wvtk == 4 )
			std::cout << "Writing tetrahedral mesh." << std::endl;
		else
		{
			std::cout << "Unknown type of mesh " << wvtk << " (1:hexahedral, 2:tetrahedral, 3:hexahedral(vtu), 4:tetrahedral(vtu)), setting hexahedral." << std::endl;
			wvtk = 1;
		}
		fvtk.open(std::string(argv[1])+"."+ext);
		if( vtu )
		{
			fvtk << "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\">" << std::endl;
			fvtk << "\t<UnstructuredGrid>" << std::endl;
		}
		else
		{
			fvtk << "# vtk DataFile Version 2.0" << std::endl;
			fvtk << "vtk file" << std::endl;
			fvtk << "ASCII" << std::endl;
			fvtk << "DATASET UNSTRUCTURED_GRID" << std::endl;
		}
	}
	
	std::cout << "intervals x " << lnx << ":" << rnx << " y " << lny << ":" << rny << " " << lnz << ":" << rnz << std::endl;
	
	std::cout << "refinement x " << refx << " y " << refy << " z " << refz << std::endl;
  
	std::cout << "Writing grid data." << std::endl;
	
	f << "DIMENS" << std::endl;
	f << (rnx-lnx)*refx << " " << (rny-lny)*refy << " " << (rnz-lnz)*refz << std::endl;
	f << "/" << std::endl;
	
	f << "SPECGRID" << std::endl;
	f << (rnx-lnx)*refx << " " << (rny-lny)*refy << " " << (rnz-lnz)*refz << " " << 1 << " " << 'P' << std::endl;
	f << "/" << std::endl;
	
	f << "COORD" << std::endl;
	for(int j = lny; j <= rny; ++j)
	{
		for(int jr = 0; jr < (j < rny ? refy : 1); jr++)
		{
			for(int i = lnx; i <= rnx; ++i)
			{
				for(int ir = 0; ir < (i < rnx ? refx : 1); ir++)
				{
					//bottom point
					xyz[0] = 240.0 * (i * 1. * refx + ir) / ( 60.0 * refx);
					xyz[1] = 440.0 * (j * 1. * refy + jr) / (220.0 * refy);
					xyz[2] = 340.0 * 0.0 / 85.0;
					transform(map,xyz,max,min,ztop,zbottom,nrmtop,nrmbottom);
					changexyz(xyz,max,min,ztop,zbottom,xyzout);
					f << xyzout[0] << " " << xyzout[1] << " " << xyzout[2];
					f << " ";
					//top point
					xyz[0] = 240.0 * (i * 1. * refx + ir) / ( 60.0 * refx);
					xyz[1] = 440.0 * (j * 1. * refy + jr) / (220.0 * refy);
					xyz[2] = 340.0 * nz / 85.0;
					transform(map,xyz,max,min,ztop,zbottom,nrmtop,nrmbottom);
					changexyz(xyz,max,min,ztop,zbottom,xyzout);
					f << xyzout[0] << " " << xyzout[1] << " " << xyzout[2];
					f << std::endl;
				}
			}
		}
	}
	f << "/" << std::endl;
	
	
	
	if( wvtk )
	{
		std::cout << "write coordinates to VTK file" << std::endl;
		size_t npx = (rnx-lnx)*refx+1;
		size_t npy = (rny-lny)*refy+1;
		size_t npz = (rnz-lnz)*refz+1;
		size_t npoints = npx*npy*npz;
		if( vtu )
		{
			size_t ncx = (rnx-lnx)*refx;
			size_t ncy = (rny-lny)*refy;
			size_t ncz = (rnz-lnz)*refz;
			size_t ncells = ncx*ncy*ncz;
			fvtk << "\t\t<Piece NumberOfPoints=\"" << npoints << "\" NumberOfCells=\"" << ncells << "\">" << std::endl;
		}
		if( vtu )
		{
			fvtk << "\t\t\t<Points>" << std::endl;
			fvtk << "\t\t\t\t<DataArray type=\"Float64\" Name=\"Points\" NumberOfComponents=\"3\" format=\"ascii\">" << std::endl;
		}
		else fvtk << "POINTS " << npoints << " double" << std::endl;
		//~ write_points_vtk(lnx,rnx,refx,lny,rnx,refy,lnz,rnz,refz,max,min);
		for(int k = lnz; k <= rnz; ++k)
		{
			for(int kr = 0; kr < (k < rnz ? refz : 1); ++kr)
			{
				for(int j = lny; j <= rny; ++j)
				{
					for(int jr = 0; jr < (j < rny ? refy : 1); jr++)
					{
						for(int i = lnx; i <= rnx; ++i)
						{
							for(int ir = 0; ir < (i < rnx ? refx : 1); ir++)
							{
								//bottom point
								xyz[0] = 240.0 * (i * 1. * refx + ir) / ( 60.0 * refx);
								xyz[1] = 440.0 * (j * 1. * refy + jr) / (220.0 * refy);
								xyz[2] = 340.0 * (k * 1. * refz + kr) / ( 85.0 * refz);
								transform(map,xyz,max,min,ztop,zbottom,nrmtop,nrmbottom);
								changexyz(xyz,max,min,ztop,zbottom,xyzout);
								fvtk << xyzout[0] << " " << xyzout[1] << " " << xyzout[2] << std::endl;
							}
						}
					}
				}
			}
		}
		fvtk << std::endl;
		if( vtu )
		{
			fvtk << "\t\t\t\t</DataArray>" << std::endl;
			fvtk << "\t\t\t</Points>" << std::endl;
		}
		std::cout << "done with coordinates in VTK file" << std::endl;
	}
	
	
	
	f << "ZCORN" << std::endl;
	nout = 0;
	for(int k = lnz; k < rnz; ++k)
	{
		for(int kr = 0; kr < refz; ++kr)
		{
			//top corners
			xyz[2] = 340.0 * (k * 1. * refz + kr) / (85.0 * refz);
			for(int j = lny; j < rny; ++j)
			{
				for(int jr = 0; jr < refy; ++jr)
				{
					xyz[1] = 440.0 * (j * 1. * refy + jr) / (220.0 * refy);
					//top corners, near left and near right
					for(int i = lnx; i < rnx; ++i)
					{
						for(int ir = 0; ir < refx; ++ir)
						{
							//top near left corner
							xyz[0] = 240.0 * (i * 1. * refx + ir) / (60.0 * refx);
							transform(map,xyz,max,min,ztop,zbottom,nrmtop,nrmbottom);
							changexyz(xyz,max,min,ztop,zbottom,xyzout);
							f << xyzout[2] << " ";
							//top near right corner
							xyz[0] = 240.0 * (i * 1. * refx + ir + 1) / (60.0 * refx);
							transform(map,xyz,max,min,ztop,zbottom,nrmtop,nrmbottom);
							changexyz(xyz,max,min,ztop,zbottom,xyzout);
							f << xyzout[2] << " ";
							nout++;
							if( nout % 5 == 0 ) 
								f << std::endl;
						}
					}
					xyz[1] = 440.0 * (j * 1. * refy + jr + 1) / (220.0 * refy);
					//top corners, far left and far right
					for(int i = lnx; i < rnx; ++i)
					{
						for(int ir = 0; ir < refx; ++ir)
						{
							// top far left corner
							xyz[0] = 240.0 * (i * 1. * refx + ir) / (60.0 * refx);
							transform(map,xyz,max,min,ztop,zbottom,nrmtop,nrmbottom);
							changexyz(xyz,max,min,ztop,zbottom,xyzout);
							f << xyzout[2] << " ";
							// top far right corner
							xyz[0] = 240.0 * (i * 1. * refx + ir + 1) / (60.0 * refx);
							transform(map,xyz,max,min,ztop,zbottom,nrmtop,nrmbottom);
							changexyz(xyz,max,min,ztop,zbottom,xyzout);
							f << xyzout[2] << " ";
							nout++;
							if( nout % 5 == 0 ) 
								f << std::endl;
						}
					}
				}
			}
			xyz[2] = 340.0 * (k * 1. * refz + kr + 1) / (85.0 * refz);
			//bottom corners 
			for(int j = lny; j < rny; ++j)
			{
				for(int jr = 0; jr < refy; ++jr)
				{
					xyz[1] = 440.0 * (j * 1. * refy + jr) / (220.0 * refy);
					//top corners, near left and near right
					for(int i = lnx; i < rnx; ++i)
					{
						for(int ir = 0; ir < refx; ++ir)
						{
							//bottom near left corner
							xyz[0] = 240.0 * (i * 1. * refx + ir) / (60.0 * refx);
							transform(map,xyz,max,min,ztop,zbottom,nrmtop,nrmbottom);
							changexyz(xyz,max,min,ztop,zbottom,xyzout);
							f << xyzout[2] << " ";
							//bottom near right corner
							xyz[0] = 240.0 * (i * 1. * refx + ir + 1) / (60.0 * refx);
							transform(map,xyz,max,min,ztop,zbottom,nrmtop,nrmbottom);
							changexyz(xyz,max,min,ztop,zbottom,xyzout);
							f << xyzout[2] << " ";
							nout++;
							if( nout % 5 == 0 ) 
								f << std::endl;
						}
					}
					xyz[1] = 440.0 * (j * 1. * refy + jr + 1) / (220.0 * refy);
					//top corners, far left and far right
					for(int i = lnx; i < rnx; ++i)
					{
						for(int ir = 0; ir < refx; ++ir)
						{
							//bottom far left corner
							xyz[0] = 240.0 * (i * 1. * refx + ir) / (60.0 * refx);
							transform(map,xyz,max,min,ztop,zbottom,nrmtop,nrmbottom);
							changexyz(xyz,max,min,ztop,zbottom,xyzout);
							f << xyzout[2] << " ";
							//bottom far right corner
							xyz[0] = 240.0 * (i * 1. * refx + ir + 1) / (60.0 * refx);
							transform(map,xyz,max,min,ztop,zbottom,nrmtop,nrmbottom);
							changexyz(xyz,max,min,ztop,zbottom,xyzout);
							f << xyzout[2] << " ";
							nout++;
							if( nout % 5 == 0 ) 
								f << std::endl;
						}
					}
				}
			}
		}
	}
	if( nout % 5 != 0 ) 
		f << std::endl;
	f << "/" << std::endl;
	
	
	if( wvtk )
	{
		std::cout << "write cells to VTK file" << std::endl;
		size_t npx = (rnx-lnx)*refx+1;
		size_t npy = (rny-lny)*refy+1;
		size_t npz = (rnz-lnz)*refz+1;
		size_t ncx = (rnx-lnx)*refx;
		size_t ncy = (rny-lny)*refy;
		size_t ncz = (rnz-lnz)*refz;
		size_t ncells = ncx*ncy*ncz;
		size_t records = 9*ncells;
		if( wvtk == 2 )
		{
			ncells *= 6;
			records = ncells*5;
		}
		if( vtu )
		{
			fvtk << "\t\t\t<Cells>" << std::endl;
			fvtk << "\t\t\t\t<DataArray type=\"UInt64\" Name=\"connectivity\" format=\"ascii\">" << std::endl;
		}
		else fvtk << "CELLS " << ncells << " " << records << std::endl;
		for(int k = lnz; k < rnz; ++k)
		{
			for(int kr = 0; kr < refz; ++kr)
			{
				for(int j = lny; j < rny; ++j)
				{
					for(int jr = 0; jr < refy; ++jr)
					{
						for(int i = lnx; i < rnx; ++i)
						{
							for(int ir = 0; ir < refx; ++ir)
							{
								size_t nvtx[8] =
								{
									(i-lnx) * refx + ir +     ((j-lny) * refy + jr    )*npx + ((k-lnz) * refz + kr    )*npx*npy,
									(i-lnx) * refx + ir + 1 + ((j-lny) * refy + jr    )*npx + ((k-lnz) * refz + kr    )*npx*npy,
									(i-lnx) * refx + ir + 1 + ((j-lny) * refy + jr + 1)*npx + ((k-lnz) * refz + kr    )*npx*npy,
									(i-lnx) * refx + ir +     ((j-lny) * refy + jr + 1)*npx + ((k-lnz) * refz + kr    )*npx*npy,
									(i-lnx) * refx + ir +     ((j-lny) * refy + jr    )*npx + ((k-lnz) * refz + kr + 1)*npx*npy,
									(i-lnx) * refx + ir + 1 + ((j-lny) * refy + jr    )*npx + ((k-lnz) * refz + kr + 1)*npx*npy,
									(i-lnx) * refx + ir + 1 + ((j-lny) * refy + jr + 1)*npx + ((k-lnz) * refz + kr + 1)*npx*npy,
									(i-lnx) * refx + ir +     ((j-lny) * refy + jr + 1)*npx + ((k-lnz) * refz + kr + 1)*npx*npy
								};
								if( wvtk == 1 || wvtk == 3)
								{
									if( !vtu ) fvtk << 8; //for vtu goes to offsets
									for(int q = 0; q < 8; ++q)
										fvtk << " " << nvtx[q];
									fvtk << std::endl;
								}
								else if( wvtk == 2 || wvtk == 4)
								{
									int ntet[6][4] =
									{
										{0,1,3,7},
										{0,1,7,5},
										{0,5,7,4},
										{1,2,3,6},
										{1,6,3,7},
										{1,6,7,5}
									};
									for(int c = 0; c < 6; ++c)
									{
										if( !vtu ) fvtk << 4; //for vtu goes to offsets
										for(int q = 0; q < 4; ++q)
											fvtk << " " << nvtx[ntet[c][q]];
										fvtk << std::endl;
									}
									
								}
							}
						}
					}
				}
			}
		}
		fvtk << std::endl;
		if( vtu )
		{
			fvtk << "\t\t\t\t</DataArray>" << std::endl;
			fvtk << "\t\t\t\t<DataArray type=\"UInt64\" Name=\"offsets\" format=\"ascii\">" << std::endl;
			size_t offset = 0;
			for(int k = lnz; k < rnz; ++k)
			{
				for(int kr = 0; kr < refz; ++kr)
				{
					for(int j = lny; j < rny; ++j)
					{
						for(int jr = 0; jr < refy; ++jr)
						{
							for(int i = lnx; i < rnx; ++i)
							{
								for(int ir = 0; ir < refx; ++ir)
								{
									if( wvtk == 1 || wvtk == 3)
									{
										offset += 8;
										fvtk << offset;
										fvtk << std::endl;
									}
									else if( wvtk == 2 || wvtk == 4)
									{
										for(int c = 0; c < 6; ++c)
										{
											offset += 4;
											fvtk << offset; 
										}
										fvtk << std::endl;
										
									}
								}
							}
						}
					}
				}
			}
			fvtk << "\t\t\t\t</DataArray>" << std::endl;
			fvtk << "\t\t\t\t<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">" << std::endl;
		}
		else fvtk << "CELL_TYPES " << ncells << std::endl;
		for(int k = lnz; k < rnz; ++k)
		{
			for(int kr = 0; kr < refz; ++kr)
			{
				for(int j = lny; j < rny; ++j)
				{
					for(int jr = 0; jr < refy; ++jr)
					{
						for(int i = lnx; i < rnx; ++i)
						{
							for(int ir = 0; ir < refx; ++ir)
							{
								int cells = 1;
								int ctype = 12;
								if( wvtk == 2 || wvtk == 4)
								{
									cells = 6;
									ctype = 10;
								}
								for(int q = 0; q < cells; ++q)
									fvtk << ctype << std::endl;
							}
						}
					}
				}
			}
		}
		fvtk << std::endl;
		if( vtu )
		{
			fvtk << "\t\t\t\t</DataArray>" << std::endl;
			fvtk << "\t\t\t</Cells>" << std::endl;
		}
		std::cout << "done with cells in VTK file" << std::endl;
	}
	
	
	std::cout << "Writing properties data." << std::endl;
	
	if( wvtk && !vtu )
	{
		size_t ncx = (rnx-lnx)*refx;
		size_t ncy = (rny-lny)*refy;
		size_t ncz = (rnz-lnz)*refz;
		size_t ncells = ncx*ncy*ncz;
		if( wvtk == 2 ) ncells *= 6;
		fvtk << "CELL_DATA " << ncells << std::endl;
	}
	
	if( wvtk && vtu )
		fvtk << "\t\t\t<CellData>" << std::endl;
	
	
	std::ifstream fporo("spe_phi.dat");
	if( fporo.fail() ) 
	{
		std::cout << "SPE10 data files not found." << std::endl;
		std::cout << "Expecting spe_phi.dat." << std::endl;
	}
	else
	{
		f << "PORO" << std::endl;
		std::vector<double> poro;
		poro.reserve(nz*ny*nx);
		for(int k = 0; k < nz; ++k)
		{
			for(int j = 0; j < ny; ++j)
			{
				for(int i = 0; i < nx; ++i)
				{
					fporo >> value;
					poro.push_back(value);
				}
			}
		}
		fporo.close();
		for(int k = 0; k < nz; ++k)
		{
			for(int kr = 0; kr < refz; ++kr)
			for(int j = 0; j < ny; ++j)
			{
				for(int jr = 0; jr < refy; ++jr)
				for(int i = 0; i < nx; ++i)
				{
					for(int ir = 0; ir < refx; ++ir)
					if( i >= lnx && i < rnx &&
						j >= lny && j < rny &&
						k >= lnz && k < rnz )
					{
						int ind = i + j * nx + k * nx*ny;
						f << poro[ind] << " ";
						nout++;
						if( nout % 10 == 0 ) 
							f << std::endl;
					}
				}
			}
		}
		if( nout % 10 != 0 ) 
			f << std::endl;
		f << "/" << std::endl;
		
		if( wvtk )
		{
			std::cout << "write PORO to VTK file" << std::endl;
			if( vtu )
			{
				fvtk << "\t\t\t\t<DataArray Name=\"PORO\" NumberOfComponents=\"1\" type=\"Float64\" format=\"ascii\">" << std::endl;
			}
			else
			{
				fvtk << "SCALARS PORO double" << std::endl;
				fvtk << "LOOKUP_TABLE default" << std::endl;
			}
			for(int k = 0; k < nz; ++k)
			{
				for(int kr = 0; kr < refz; ++kr)
				for(int j = 0; j < ny; ++j)
				{
					for(int jr = 0; jr < refy; ++jr)
					for(int i = 0; i < nx; ++i)
					{
						for(int ir = 0; ir < refx; ++ir)
						if( i >= lnx && i < rnx &&
							j >= lny && j < rny &&
							k >= lnz && k < rnz )
						{
							int ind = i + j * nx + k * nx*ny;
							int cells = 1;
							if( wvtk == 2 ) cells = 6;
							for(int q = 0; q < cells; ++q)
								fvtk << poro[ind] << std::endl;
						}
					}
				}
			}
			fvtk << std::endl;
			if( vtu )
				fvtk << "\t\t\t\t</DataArray>" << std::endl;
			std::cout << "done with PORO in VTK file" << std::endl;
		}
	}
	
	
	std::ifstream fperm("spe_perm.dat");
	if( fperm.fail() ) 
	{
		std::cout << "SPE10 data files not found." << std::endl;
		std::cout << "Expecting spe_perm.dat." << std::endl;
		return false;
	}
	else
	{
		std::vector<double> perm[3], permnew[6];
		for(int l = 0; l < 3; ++l)
		{
			perm[l].reserve(nz*ny*nx);
			for(int k = 0; k < nz; ++k)
			{
				for(int j = 0; j < ny; ++j)
				{
					for(int i = 0; i < nx; ++i)
					{
						fperm >> value;
						perm[l].push_back(value);
					}
				}
			}
		}
		fperm.close();
		
		for(int l = 0; l < 6; ++l)
			permnew[l].reserve(nz*ny*nz);
		
		nout = 0;
		for(int k = 0; k < nz; ++k)
		{
			for(int j = 0; j < ny; ++j)
			{
				for(int i = 0; i < nx; ++i)
				{
					if( i >= lnx && i < rnx &&
						j >= lny && j < rny &&
						k >= lnz && k < rnz )
					{
						xyz[0] = 240.0 * (i+0.5) / 60.0;
						xyz[1] = 440.0 * (j+0.5) / 220.0;
						xyz[2] = 340.0 * (k+0.5) / 85.0;
						Kin[0] = perm[0][nout];
						Kin[1] = perm[1][nout];
						Kin[2] = perm[2][nout];
						transform(map,xyz,max,min,ztop,zbottom,nrmtop,nrmbottom);
						for(int l = 0; l < 3; ++l)
							nrm[l] = (nrmtop[l]-nrmbottom[l])*(xyz[2]-min[2])/(max[2]-min[2]) + nrmbottom[l];
						rotate_tensor(nrm,Kin,Kout);
						for(int l = 0; l < 6; ++l)
							permnew[l].push_back(Kout[l]);
					}
					else for(int l = 0; l < 6; ++l)
						permnew[l].push_back(0.0);
					nout++;
				}
			}
		}
		
		f << "MPFA" << std::endl;
		f << 1 << " " << 0 << std::endl; // define tensor in x,y,z coords and use TPFA
		f << "/" << std::endl;
		const char c1[6] = {'X','X','X','Y','Y','Z'};
		const char c2[6] = {'X','Y','Z','Y','Z','Z'};
		for(int l = 0; l < 6; ++l)
		{
			f << "PERM" << c1[l] << c2[l] << std::endl;
			m = nout = 0;
			for(int k = 0; k < nz; ++k)
			{
				for(int kr = 0; kr < refz; ++kr)
				for(int j = 0; j < ny; ++j)
				{
					for(int jr = 0; jr < refy; ++jr)
					for(int i = 0; i < nx; ++i)
					{
						for(int ir = 0; ir < refx; ++ir)
						if( i >= lnx && i < rnx &&
							j >= lny && j < rny &&
							k >= lnz && k < rnz )
						{
							int ind = i + j * nx + k * nx*ny;
							f << permnew[l][ind] << " ";
							nout++;
							if( nout % 10 == 0 ) 
								f << std::endl;
						}
					}
				}
			}
			if( nout % 10 != 0 ) 
				f << std::endl;
			f << "/" << std::endl;
		}
		
		
		if( wvtk )
		{
			std::cout << "write PERM to VTK file" << std::endl;
			if( vtu )
			{
				fvtk << "\t\t\t\t<DataArray Name=\"PERM\" NumberOfComponents=\"6\" type=\"Float64\" format=\"ascii\">" << std::endl;
			}
			else
			{
				fvtk << "SCALARS PERM double 6" << std::endl;
				fvtk << "LOOKUP_TABLE default" << std::endl;
			}
			for(int k = 0; k < nz; ++k)
			{
				for(int kr = 0; kr < refz; ++kr)
				for(int j = 0; j < ny; ++j)
				{
					for(int jr = 0; jr < refy; ++jr)
					for(int i = 0; i < nx; ++i)
					{
						for(int ir = 0; ir < refx; ++ir)
						if( i >= lnx && i < rnx &&
							j >= lny && j < rny &&
							k >= lnz && k < rnz )
						{
							int ind = i + j * nx + k * nx*ny;
							int cells = 1;
							if( wvtk == 2 ) cells = 6;
							for(int q = 0; q < cells; ++q)
							{
								for(int l = 0; l < 6; ++l)
									fvtk << permnew[l][ind] << " ";
								fvtk << std::endl;
							}
						}
					}
				}
			}
			fvtk << std::endl;
			if( vtu )
				fvtk << "\t\t\t\t</DataArray>" << std::endl;
			std::cout << "done with PERM in VTK file" << std::endl;
		}
	}
	
	std::cout << "Closing file " << argv[1] << "." << std::endl;
	f.close();
	if( wvtk )
	{
		if( vtu )
		{
			fvtk << "\t\t\t</CellData>" << std::endl;
			fvtk << "\t\t</Piece>" << std::endl;
			fvtk << "\t</UnstructuredGrid>" << std::endl;
			fvtk << "</VTKFile>" << std::endl;
		}
		std::cout << "Closing VTK file" << std::endl;
		fvtk.close();
	}
	std::cout << "Done!" << std::endl;

	return 0;
}
