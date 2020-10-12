#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>


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


void transform(double xyz[3], double max[3], double min[3], double & ztop, double & zbottom, double nrmtop[3], double nrmbottom[3])
{
	double x = (xyz[0]-min[0])/(max[0]-min[0]), y = (xyz[1]-min[1])/(max[1]-min[1]), z = (xyz[2]-min[2])/(max[2]-min[2]);
	if( x < 0 || x > 1 ) std::cout << "x: " << x << std::endl;
	if( y < 0 || y > 1 ) std::cout << "y: " << y << std::endl;
	if( z < 0 || z > 1 ) std::cout << "z: " << z << " xyz " << xyz[0] << " " << xyz[1] << " " << xyz[2] << " unit " << x << " " << y << " " << z << " min " << min[0] << " " << min[1] << " " << min[2] << " max " << max[0] << " " << max[1] << " " << max[2] << std::endl;
	double dztopdx = 0, dztopdy = 0;
	double dzbottomdx = 0, dzbottomdy = 0;
	//zbottom = x*x*4-y*y*4-sin(6*x)*8 - cos(4*y)*4 - x*15;
	//ztop = x*x*4-y*y*4-sin(6*x)*8 - cos(4*y)*4 - y*15 + 15;
	zbottom = std::sin(y*pi *2)*0.2;// ((x - 0.5)*(x - 0.5) + (y - 0.5)*(y - 0.5))*0.4;
	ztop = std::sin(y*pi *2)*0.2+1;// 1 + cos(4 * y)*0.2 + ((x - 0.5)*(x - 0.5) + (y - 0.5)*(y - 0.5))*0.4 + 1 + cos(4 * y)*0.5;
	
	dzbottomdx = 0;
	dzbottomdy = 2*pi*std::cos(y*pi*2)*0.2;
	
	dztopdx = 0;
	dztopdy = 2*pi*std::cos(y*pi*2)*0.2;
	
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

int main(int argc, char *argv[]) 
{
	if (argc < 2)
	{
		std::cout << "Usage: " << argv[0] ;
		std::cout << " output.grdecl [lnx=0 rnx=60 lny=0 rny=220 lnz=0 rnz=85]" << std::endl;
		return -1;
	}

  
	std::cout << "Opening " << argv[1] << " for output." << std::endl;
	std::ofstream f(argv[1]);
	
	if( f.fail() )
	{
		std::cout << "Cannot open " << argv[1] << " for writing!" << std::endl;
		return -1;
	}
	
	double value, ztop, zbottom, xyz[3], xyzout[3], nrmtop[3], nrmbottom[3], nrm[3], Kin[3], Kout[6];
	double max[3] = {240,440,340}, min[3] = {0,0,0};
	int nx = 60, ny = 220, nz = 85, nout = 0, m;
	int lnx = 0, lny = 0, lnz = 0;
	int rnx = nx, rny = ny, rnz = nz;
	
	if( argc > 2 ) lnx = atoi(argv[2]);
	if( argc > 3 ) rnx = atoi(argv[3]);
	if( argc > 4 ) lny = atoi(argv[4]);
	if( argc > 5 ) rny = atoi(argv[5]);
	if( argc > 6 ) lnz = atoi(argv[6]);
	if( argc > 7 ) rnz = atoi(argv[7]);
	std::cout << "intervals x " << lnx << ":" << rnx << " y " << lny << ":" << rny << " " << lnz << ":" << rnz << std::endl;
  
	std::cout << "Writing grid data." << std::endl;
	
	f << "DIMENS" << std::endl;
	f << rnx-lnx << " " << rny-lny << " " << rnz-lnz << std::endl;
	f << "/" << std::endl;
	
	f << "SPECGRID" << std::endl;
	f << rnx-lnx << " " << rny-lny << " " << rnz-lnz << " " << 1 << " " << 'P' << std::endl;
	f << "/" << std::endl;
	
	f << "COORD" << std::endl;
	for(int j = lny; j <= rny; ++j)
	{
		for(int i = lnx; i <= rnx; ++i)
		{
			//bottom point
			xyz[0] = 240.0 * i / 60.0;
			xyz[1] = 440.0 * j / 220.0;
			xyz[2] = 340.0 * 0.0 / 85.0;
			transform(xyz,max,min,ztop,zbottom,nrmtop,nrmbottom);
			changexyz(xyz,max,min,ztop,zbottom,xyzout);
			f << xyzout[0] << " " << xyzout[1] << " " << xyzout[2];
			f << " ";
			//top point
			xyz[0] = 240.0 * i / 60.0;
			xyz[1] = 440.0 * j / 220.0;
			xyz[2] = 340.0 * nz / 85.0;
			transform(xyz,max,min,ztop,zbottom,nrmtop,nrmbottom);
			changexyz(xyz,max,min,ztop,zbottom,xyzout);
			f << xyzout[0] << " " << xyzout[1] << " " << xyzout[2];
			f << std::endl;
		}
	}
	f << "/" << std::endl;
	
	f << "ZCORN" << std::endl;
	nout = 0;
	for(int k = lnz; k < rnz; ++k)
	{
		//top corners
		xyz[2] = 340.0 * k / 85.0;
		for(int j = lny; j < rny; ++j)
		{
			xyz[1] = 440.0 * j / 220.0;
			//top corners, near left and near right
			for(int i = lnx; i < rnx; ++i)
			{
				//top near left corner
				xyz[0] = 240.0 * i / 60.0;
				transform(xyz,max,min,ztop,zbottom,nrmtop,nrmbottom);
				changexyz(xyz,max,min,ztop,zbottom,xyzout);
				f << xyzout[2] << " ";
				//top near right corner
				xyz[0] = 240.0 * (i+1) / 60.0;
				transform(xyz,max,min,ztop,zbottom,nrmtop,nrmbottom);
				changexyz(xyz,max,min,ztop,zbottom,xyzout);
				f << xyzout[2] << " ";
				nout++;
				if( nout % 5 == 0 ) 
					f << std::endl;
			}
			xyz[1] = 440.0 * (j+1) / 220.0;
			//top corners, far left and far right
			for(int i = lnx; i < rnx; ++i)
			{
				// top far left corner
				xyz[0] = 240.0 * i / 60.0;
				transform(xyz,max,min,ztop,zbottom,nrmtop,nrmbottom);
				changexyz(xyz,max,min,ztop,zbottom,xyzout);
				f << xyzout[2] << " ";
				// top far right corner
				xyz[0] = 240.0 * (i+1) / 60.0;
				transform(xyz,max,min,ztop,zbottom,nrmtop,nrmbottom);
				changexyz(xyz,max,min,ztop,zbottom,xyzout);
				f << xyzout[2] << " ";
				nout++;
				if( nout % 5 == 0 ) 
					f << std::endl;
			}
		}
		xyz[2] = 340.0 * (k+1) / 85.0;
		//bottom corners 
		for(int j = lny; j < rny; ++j)
		{
			xyz[1] = 440.0 * j / 220.0;
			//top corners, near left and near right
			for(int i = lnx; i < rnx; ++i)
			{
				//bottom near left corner
				xyz[0] = 240.0 * i / 60.0;
				transform(xyz,max,min,ztop,zbottom,nrmtop,nrmbottom);
				changexyz(xyz,max,min,ztop,zbottom,xyzout);
				f << xyzout[2] << " ";
				//bottom near right corner
				xyz[0] = 240.0 * (i+1) / 60.0;
				transform(xyz,max,min,ztop,zbottom,nrmtop,nrmbottom);
				changexyz(xyz,max,min,ztop,zbottom,xyzout);
				f << xyzout[2] << " ";
				nout++;
				if( nout % 5 == 0 ) 
					f << std::endl;
			}
			xyz[1] = 440.0 * (j+1) / 220.0;
			//top corners, far left and far right
			for(int i = lnx; i < rnx; ++i)
			{
				//bottom far left corner
				xyz[0] = 240.0 * i / 60.0;
				transform(xyz,max,min,ztop,zbottom,nrmtop,nrmbottom);
				changexyz(xyz,max,min,ztop,zbottom,xyzout);
				f << xyzout[2] << " ";
				//bottom far right corner
				xyz[0] = 240.0 * (i+1) / 60.0;
				transform(xyz,max,min,ztop,zbottom,nrmtop,nrmbottom);
				changexyz(xyz,max,min,ztop,zbottom,xyzout);
				f << xyzout[2] << " ";
				nout++;
				if( nout % 5 == 0 ) 
					f << std::endl;
			}
		}
	}
	if( nout % 5 != 0 ) 
		f << std::endl;
	f << "/" << std::endl;
	
	
	std::cout << "Writing properties data." << std::endl;
	
	
	std::ifstream fporo("spe_phi.dat");
	if( fporo.fail() ) 
	{
		std::cout << "SPE10 data files not found." << std::endl;
		std::cout << "Expecting spe_phi.dat." << std::endl;
	}
	else
	{
		f << "PORO" << std::endl;
		for(int k = 0; k < nz; ++k)
		{
			for(int j = 0; j < ny; ++j)
			{
				for(int i = 0; i < nx; ++i)
				{
					fporo >> value;
					if( i >= lnx && i < rnx &&
						j >= lny && j < rny &&
						k >= lnz && k < rnz )
					{
						f << value << " ";
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
		fporo.close();
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
						transform(xyz,max,min,ztop,zbottom,nrmtop,nrmbottom);
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
				for(int j = 0; j < ny; ++j)
				{
					for(int i = 0; i < nx; ++i)
					{
						if( i >= lnx && i < rnx &&
							j >= lny && j < rny &&
							k >= lnz && k < rnz )
						{
							f << permnew[l][m] << " ";
							nout++;
							if( nout % 10 == 0 ) 
								f << std::endl;
						}
						m++;
					}
				}
			}
			if( nout % 10 != 0 ) 
				f << std::endl;
			f << "/" << std::endl;
		}
	}
	
	std::cout << "Closing file " << argv[1] << "." << std::endl;
	f.close();
	std::cout << "Done!" << std::endl;

	return 0;
}
