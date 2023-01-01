/////////////////////////////////////////////////
/// @file Outputters.cpp
/// @brief Generation of output writing to disk
/////////////////////////////////////////////////

#include "Outputters.h"

namespace fs = std::filesystem;
namespace chrono = std::chrono;

/// Directly generate the simCode at construction
OutputManager::OutputManager() {
  simCode = SimCodeGenerator();
  outputStyle = 'c';
}

/// Generate the identifier number reverse from year to minute in the format
/// yy-mm-dd_hh-MM-ss
std::string OutputManager::SimCodeGenerator() {
  const chrono::time_point<chrono::system_clock> now{
      chrono::system_clock::now()};
  const chrono::year_month_day ymd{chrono::floor<chrono::days>(now)};
  const auto tod = now - chrono::floor<chrono::days>(now);
  const chrono::hh_mm_ss hms{tod};

  std::stringstream temp;
  temp << std::setfill('0') << std::setw(2)
       << static_cast<int>(ymd.year() - chrono::years(2000)) << "-"
       << std::setfill('0') << std::setw(2)
       << static_cast<unsigned>(ymd.month()) << "-"
       << std::setfill('0') << std::setw(2)
       << static_cast<unsigned>(ymd.day()) << "_"
       << std::setfill('0') << std::setw(2) << hms.hours().count()
       << "-" << std::setfill('0')
       << std::setw(2) << hms.minutes().count() << "-"
       << std::setfill('0') << std::setw(2)
       << hms.seconds().count();
       //<< "_" << hms.subseconds().count(); // subseconds render the filename
       // too large
  return temp.str();
}

/** Generate the folder to save the data to by one process:
 * In the given directory it creates a direcory "SimResults" and a directory
 * with the simCode. The relevant part of the main file is written to a
 * "config.txt" file in that directory to log the settings. */
void OutputManager::generateOutputFolder(const std::string &dir) {
  // Do this only once for the first process
  int myPrc = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &myPrc);
  if (myPrc == 0) {
    if (!fs::is_directory(dir))
      fs::create_directory(dir);
    if (!fs::is_directory(dir + "/SimResults"))
      fs::create_directory(dir + "/SimResults");
    if (!fs::is_directory(dir + "/SimResults/" + simCode))
      fs::create_directory(dir + "/SimResults/" + simCode);
  }
  // path variable for the output generation
  Path = dir + "/SimResults/" + simCode + "/";

  // Logging configurations from main.cpp -> no more necessary
  /*
  std::ifstream fin("main.cpp");
  std::ofstream fout(Path + "config.txt");
  std::string line;
  int begin = 1000;
  for (int i = 1; !fin.eof(); i++) {
    getline(fin, line);
    if (line.starts_with("    //------------ B")) {
        begin=i;
    }
    if (i < begin) {
      continue;
    }
    fout << line << std::endl;
    if (line.starts_with("    //------------- E")) {
        break;
    }
  }
  */
  return;
}

void OutputManager::set_outputStyle(const char _outputStyle){
    outputStyle = _outputStyle;
}

/** Write the field data either in csv format to one file per each process
 * (patch) or in binary form to a single file. Files are stores inthe simCode
 * directory. For csv files the state (simulation step) denotes the
 * prefix and the suffix after an underscore is given by the process/patch
 * number. Binary files are simply named after the step number. */
void OutputManager::outUState(const int &state, const Lattice &lattice,
        const LatticePatch &latticePatch) {
  switch(outputStyle){
      case 'c': { // one csv file per process
                    std::ofstream ofs;
  ofs.open(Path + std::to_string(state) + "_"
          + std::to_string(lattice.my_prc) + ".csv");
  // Precision of sunrealtype in significant decimal digits; 15 for IEEE double
  ofs << std::setprecision(std::numeric_limits<sunrealtype>::digits10);

  // Walk through each lattice point
  const sunindextype totalNP = latticePatch.discreteSize();
  for (sunindextype i = 0; i < totalNP * 6; i += 6) {
    // Six columns to contain the field data: Ex,Ey,Ez,Bx,By,Bz
    ofs << latticePatch.uData[i + 0] << "," << latticePatch.uData[i + 1] << ","
        << latticePatch.uData[i + 2] << "," << latticePatch.uData[i + 3] << ","
        << latticePatch.uData[i + 4] << "," << latticePatch.uData[i + 5]
        << std::endl;
  }
  ofs.close();
  break;
                }

      case 'b': { // a single binary file
  // Open the output file
  MPI_File fh;
  const std::string filename = Path+std::to_string(state);
  MPI_File_open(lattice.comm,&filename[0],MPI_MODE_WRONLY|MPI_MODE_CREATE,
          MPI_INFO_NULL,&fh);
  // number of datapoints in the patch with process offset
  const sunindextype count = latticePatch.discreteSize()*
      lattice.get_dataPointDimension();
  MPI_Offset offset = lattice.my_prc*count*sizeof(MPI_SUNREALTYPE);
  // Go to offset in file and write data to it; maximal precision in
  // "native" representation
  MPI_File_set_view(fh,offset,MPI_SUNREALTYPE,MPI_SUNREALTYPE,"native",
          MPI_INFO_NULL);
  MPI_Request write_request;
  MPI_File_iwrite_all(fh,latticePatch.uData,count,MPI_SUNREALTYPE,
          &write_request);
  MPI_Wait(&write_request,MPI_STATUS_IGNORE);
  MPI_File_close(&fh);
  break;
                }
      default: {
  errorKill("No valid output style defined."
          " Choose between (c): one csv file per process,"
          " (b) one binary file");
  break;
               }
  }
}

