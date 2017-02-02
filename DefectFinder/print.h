// How to use: Include this file in your main.cpp (#include "path/to/file/print.h").
// In this file: Edit the constant string in line 35 to your output folder
// In main:
// Create an object Print print; and call the initializer print.init(ID, noCells)
//   The arguments are:
//   - ID: The name of the simulation as a string e.g. "0042"
//   - noCells: The number of particles
//
//   It will create the folder to write to in 'path' + 2 subfolders for normal data
//   files and the visualisation. It also copies vid.html to path+run+"/vid/vid.html"
//   Finally call print_data_js(x1, y1, x2, y2, r [, c])   [last one is optional]
//   This will print the x, y position of both disks of the bacterium, the radius of
//   the disk, and c to change the colour of the cell. Colour should be an integer
//   between 0 and 360. This number corresponds to the hue colour
//
// To write cell positions to a file, call print.print_data_js(x1, y1, x2, y2, r)
// or print_data_js(x1, y1, x2, y2, r, colour)

#ifdef _WIN32
    const int unixSystem=0;
#elif __linux__
    const int unixSystem=1;
#else
    const int unixSystem=1;
#endif

#include <iostream> // For cout
#include <iomanip>  // For setprecision
#include <fstream>  // For ofstream
#include <string>   // For string
#include <cstdlib>  // For exit()

using namespace std;

const string path = "C:\\Users\\Duco\\Documents\\1Studie\\1MEP\\Code\\DefectFinder\\Visual\\";
// On windows this will look like
// path = "C:\\Users\\name\\output folder\\";
// It doesn't matter if your path contains spaces
// On linux this will look like
// path = "/home/name/outputFolder/";

struct Print
{
    Print();
    ~Print();

    void init(string, int);                 // Initialize filesystem and output streams
    void set_no_particles(int);             // Use this to set the number of particles if it changes
	void print_data_js(double, double, double, double, double, int, int);
    void print_data_js(double, double, double, double, double, int);
    void print_data_js(double, double, double, double, double);
    void print_data(long int, double, double);
    void resize(double &, double &, double &, double &, double &);

    bool first;
    int N;
    int p;
    string run;
    ofstream data_js, data;
};

Print::Print()
{
    // Empty
}

Print::~Print()
{
	data_js << "]\n];";

    data_js.close();
    data.close();
}

void Print::init(string ID, int noCells)
// Sets up the structure in the filesystem to save the data and opens the files for writing
// path /
//        ID /
//             dat /
//                   data.dat
//             vid /
//                   vid.html
//                   data.js
//        jquery.js (make sure to put that file in the folder specified by path)      
{
    first = true;
    N = noCells;
    p = 0;
    run = ID;
    int check = 0;

    if(unixSystem == 1) // Use '/' for filesystem, doesn't need quotes for files/folders with spaces, use cp to copy
    {
        ifstream test;
        test.open((path+run+"/vid/vid.html").c_str());  // Check if this is gonna overwrite data

        if(test.fail())
        {
            check += system(("mkdir "+path+run).c_str());
            check += system(("mkdir "+path+run+"/dat").c_str());
            check += system(("mkdir "+path+run+"/vid").c_str());
            if(check != 0)
            {
                cout << "Failed to create folders. Status 715\n"; exit(715);
            }

            check = 0;

            check += system(("cp "+path+"vid.html "+path+run+"/vid/vid.html").c_str());

            if(check != 0)
            {
                cout << "Failed to copy files to vid/. Status 716\n"; exit(716);
            }
        }
        else
        {
            test.close();
            cout << "This action will overwrite old data. Press \'9\' to continue.\n";
            cin >> check;
            if(check!=9)
                exit(9);
        }

        data_js.open((path+run+"/vid/data.js").c_str());
        if(data_js.fail())
        {
            cout << "Failed to open " << path+run << "/vid/data.js. Status 717\n";
            exit(717);
        }
        data.open((path+run+"/dat/data.dat").c_str());
        if(data.fail())
        {
            cout << "Failed to open " << path+run << "/vid/data.dat. Status 718\n";
            exit(718);
        }
    }
    else        // Use '\\' for filesystem, needs quotes for files/folders with spaces, use COPY to copy
    {
        ifstream test;
        test.open((path+run+"\\vid\\vid.html").c_str());  // Check if this is gonna overwrite data

        if(test.fail())
        {   // Don't actually do these tests because windows seems to screw this up
            check += system(("mkdir \""+path+run+"\"").c_str());
            check += system(("mkdir \""+path+run+"\\dat\"").c_str());
            check += system(("mkdir \""+path+run+"\\vid\"").c_str());
            check += system(("copy \""+path+"vid.html\" \""+path+run+"\\vid\\vid.html\"").c_str());
            if(check != 0)
            {
                cout << "Failed to create folders. Status 715\n"; exit(715);
            }

            check = 0;

            if(check != 0)
            {
                cout << "Failed to copy files to vid/. Status 716\n"; exit(716);
            }
        }
        else
        {
            test.close();
            cout << "This action will overwrite old data. Press \'9\' to continue.\n";
            cin >> check;
            if(check!=9)
                exit(9);
        }

        data_js.open((path+run+"\\vid\\data.js").c_str());
        if(data_js.fail())
        {
            cout << "Failed to open " << (path+run+"\\vid\\data.js").c_str() << ". Status 717\n";
            exit(717);
        }
        data.open((path+run+"\\dat\\data.dat").c_str());
        if(data.fail())
        {
            cout << "Failed to open " << (path+run+"\\vid\\data.dat").c_str() << ". Status 718\n";
            exit(718);
        }
    }

    data_js << "var cells = [\n";
}

void Print::set_no_particles(int noParticles)
// Set the number of particles (Alternatively you can use set N with print.N = noParticles)
{
    N = noParticles;
}

void Print::print_data_js(double x1, double y1, double x2, double y2, double r, int c, int id)
// Prints the information about the cell to the data file
// When called, it prints the x- and y-position of both disks, the radius of the disk, 
// an integer representing the color (could indicate overlap) and an integer containing the id of the cell (cell nr)
{
	++p;

	if (p == 1 && first)
	{
		data_js << "[";
		first = false;
	}
	else if (p == 1)
		data_js << "],\n[";

	resize(x1, y1, x2, y2, r);
	int ix1 = x1;
	int iy1 = y1;
	int ix2 = x2;
	int iy2 = y2;

	data_js << "[" << ix1 << "," << iy1 << "," << ix2 << "," << iy2 << "," << r << "," << c << "," << id << "]";

	if (p == N)
	{
		data_js << "\n";
		p = 0;
	}
	else
		data_js << ",\n";
}

void Print::print_data_js(double x1, double y1, double x2, double y2, double r, int c)
// Prints the information about the cell to the data file
// When called, it prints the x- and y-position of both disks, the radius of the disk, 
// and an integer representing the color (could indicate overlap)
{
    ++p;

    if(p==1 && first)
    {
        data_js << "[";
        first = false;
	}
    else if(p==1)
        data_js << "],\n[";

    resize(x1, y1, x2, y2, r);
    int ix1 = x1;
    int iy1 = y1;
    int ix2 = x2;
    int iy2 = y2;

    data_js << "[" << ix1 << "," << iy1 << "," << ix2 << "," << iy2 << "," << r << "," << c << "]";

    if(p==N)
    {
        data_js << "\n";
        p=0;
    }
    else
        data_js << ",\n";
}

void Print::print_data_js(double x1, double y1, double x2, double y2, double r)
// Prints the information about the cell to the data file
// When called, it prints the x- and y-position of both disks, the radius of the disk, 
// and an integer representing the color (could indicate overlap)
{
    ++p;

    if(p==1 && first)
    {
        data_js << "[";
        first = false;
	}
    else if(p==1)
        data_js << "],\n[";

    resize(x1, y1, x2, y2, r);
    int ix1 = x1;
    int iy1 = y1;
    int ix2 = x2;
    int iy2 = y2;

    data_js << "[" << ix1 << "," << iy1 << "," << ix2 << "," << iy2 << "," << r << "]";

    if(p==N)
    {
        data_js << "\n";
        p=0;
    }
    else
        data_js << ",\n";
}

void Print::print_data(long int time, double firstQuantity, double secondQuantity)
// Print data to tab separated datafile every XXX simulation steps. Useful for plotting data in time
// You can add or remove quantities that you want to make plots of.
{
    data << time << "\t" << firstQuantity << "\t" << secondQuantity << "\n";
}

void Print::resize(double &x1, double &y1, double &x2, double &y2, double &r)
// Resizes (and translates) the cells such that they are beautifully displayed
{
    x1 *= 20;
    y1 *= -20;
    x1 += 1000;
    y1 += 1000;
    x2 *= 20;
    y2 *= -20;
    x2 += 1000;
    y2 += 1000;
    r *= 20;
}
