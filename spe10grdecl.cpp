#include <iostream>
#include <fstream>

int main(int argc, char *argv[]) 
{
	if (argc < 2)
	{
		std::cout << "Usage: " << argv[0] ;
		std::cout << " [output=spe10.grdecl]" << std::endl;
		return -1;
	}

  
	std::ofstream f(argv[1]);
	double value;
	int nx = 60, ny = 220, nz = 85, nout = 0;
  
	
	f << "DIMENS" << std::endl;
	f << nx << " " << ny << " " << nz << std::endl;
	f << "/" << std::endl;
	
	f << "SPECGRID" << std::endl;
	f << nx << " " << ny << " " << nz << " " << 1 << " " << 'P' << std::endl;
	f << "/" << std::endl;
	
	f << "COORD" << std::endl;
	for(int i = 0; i <= nx; ++i)
	{
		for(int j = 0; j <= ny; ++j)
		{
			//bottom point
			f << 240.0 * i / 60.0 << " " << 240.0 * j / 60.0 << " " << 340.0 * 0.0 / 85.0;
			f << " ";
			//top point
			f << 240.0 * i / 60.0 << " " << 240.0 * j / 60.0 << " " << 340.0 * 1.0 / 85.0;
			f << std::endl;
		}
	}
	f << "/" << std::endl;
	
	f << "ZCORN" << std::endl;
	for(int i = 0; i < nx; ++i)
	{
		for(int j = 0; j < ny; ++j)
		{
			for(int k = 0; k < nz; ++k)
			{
				for(int l = 0; l < 8; ++l)
				{
					
				}
			}
		}
	}
	f << "/" << std::endl;
	
	
	std::ifstream fporo("spe_phi.dat");
	if( poro.fail() ) 
	{
		std::cout << "SPE10 data files not found." << std::endl;
		std::cout << "Expecting spe_perm.dat and spe_phi.dat." << std::endl;
		return false;
	}
	f << "PORO" << std::endl;
	for(int k = 0; k < nz; ++k)
	{
		for(int j = 0; j < ny; ++j)
		{
			for(int i = 0; i < nx; ++i)
			{
				fporo >> value;
				f << value << " ";
				nout++;
				if( nout % 10 ) 
					f << std::endl;
			}
		}
	}
	f << "/" << std::endl;
	fporo.close();
	
	
	std::ifstream fperm("spe_perm.dat");
	if( fperm.fail() ) 
	{
		std::cout << "SPE10 data files not found." << std::endl;
		std::cout << "Expecting spe_perm.dat and spe_phi.dat." << std::endl;
		return false;
	}
	const char c[3] = {'X','Y','Z'};
	for(int l = 0; l < 3; ++l)
	{
		f << "PERM" << c[l] << std::endl;
		nout = 0;
		for(int k = 0; k < nz; ++k)
		{
			for(int j = 0; j < ny; ++j)
			{
				for(int i = 0; i < nx; ++i)
				{
					fperm >> value;
					f << value << " ";
					nout++;
					if( nout % 10 ) 
						f << std::endl;
				}
			}
		}
		f << "/" << std::endl;
	}
	fperm.close();
	
	f.close();

	return 0;
}
