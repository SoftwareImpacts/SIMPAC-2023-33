//////////////////////////////////////////////////////////////$
/// @file Outputters.cpp
/// @brief Generation of output writing to disk
//////////////////////////////////////////////////////////////$

#include "Outputters.h"

/// Directly generate the simCode at construction
OutputManager::OutputManager() {
  simCode = SimCodeGenerator();
  MPI_Comm_rank(MPI_COMM_WORLD, &myPrc);
}

/// Generate the identifier number reverse from year to minute in the format
/// yy-mm-dd_hh-MM-ss
string OutputManager::SimCodeGenerator() {
  const chrono::time_point<chrono::system_clock> now{
      chrono::system_clock::now()};
  const chrono::year_month_day ymd{chrono::floor<chrono::days>(now)};
  const auto tod = now - chrono::floor<chrono::days>(now);
  const chrono::hh_mm_ss hms{tod};

  stringstream temp;
  temp << setfill('0') << setw(2)
       << static_cast<int>(ymd.year() - chrono::years(2000)) << "-"
       << setfill('0') << setw(2) << static_cast<unsigned>(ymd.month()) << "-"
       << setfill('0') << setw(2) << static_cast<unsigned>(ymd.day()) << "_"
       << setfill('0') << setw(2) << hms.hours().count() << "-" << setfill('0')
       << setw(2) << hms.minutes().count() << "-" << setfill('0') << setw(2)
       << hms.seconds().count();
  //<< "_" << hms.subseconds().count(); // subseconds render the filename too
  //large
  return temp.str();
}

/** Generate the folder to save the data to by one process:
 * In the given directory it creates a direcory "SimResults" and a directory
 * with the simCode. The relevant part of the main file is written to a
 * "config.txt" file in that directory to log the settings. */
void OutputManager::generateOutputFolder(const string &dir) {
  // Do this only once for the first process
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

  // Logging configurations from main.cpp
  ifstream fin("main.cpp");
  ofstream fout(Path + "config.txt");
  string line;
  int begin=1000;
  for (int i = 1; !fin.eof(); i++) {
    getline(fin, line);
    if (line.starts_with("    //------------ B")) {
        begin=i;
    }
    if (i < begin) {
      continue;
    }
    fout << line << endl;
    if (line.starts_with("    //------------- E")) {
        break;
    }
  }

  return;
}

/** Write the field data to a csv file from each process (patch) with the field
 * data into the simCode directory. The state (simulation step) denotes the
 * prefix and the suffix after an underscore is given by the process/patch
 * number */
void OutputManager::outUState(const int &state, const LatticePatch &latticePatch) {
  ofstream ofs;
  ofs.open(Path + to_string(state) + "_" + to_string(myPrc) + ".csv");
  // Set precision, number of digits for the values
  ofs << setprecision(numeric_limits<sunrealtype>::digits10);

  // Walk through each lattice point
  for (int i = 0; i < latticePatch.discreteSize() * 6; i += 6) {
    // Six columns to contain the field data: Ex,Ey,Ez,Bx,By,Bz
    ofs << latticePatch.uData[i + 0] << "," << latticePatch.uData[i + 1] << ","
        << latticePatch.uData[i + 2] << "," << latticePatch.uData[i + 3] << ","
        << latticePatch.uData[i + 4] << "," << latticePatch.uData[i + 5]
        << endl;
  }

  ofs.close();

  return;
}

/// Return the date+time simulation identifier for logging
string OutputManager::getSimCode() { return simCode; }
