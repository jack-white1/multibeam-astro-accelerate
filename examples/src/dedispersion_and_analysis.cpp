/**
 * Example code for linking against astro-accelerate library.
 *
 * Compile with: g++ -std=c++11 -I/path/to/astro-accelerate/include/ -L/path/to/astro-accelerate/build -Wl,-rpath,/path/to/astro-accelerate/build -lastroaccelerate -I/usr/local/cuda/include/ -L/usr/local/cuda/lib64 -Wl,-rpath,/usr/local/cuda/lib64 -I/usr/local/cuda-8.0/samples/common/inc/ -lcudart dedispersion.cpp -o test
 */

#include <sys/shm.h>
#include <iostream>
#include <fstream>
#include <sys/time.h>
#include <bits/stdc++.h>
#include <ctime>
#include <sys/stat.h> // new edits, to create a directory; Added by Rushikesh
#include "aa_ddtr_plan.hpp"
#include "aa_ddtr_strategy.hpp"
#include "aa_filterbank_metadata.hpp"
#include "aa_sigproc_input.hpp"
#include "aa_permitted_pipelines_generic.hpp"
#include "aa_pipeline_api.hpp"
#include "aa_device_info.hpp"
#include "gmrt_newcorr.h"
#include "frb_shm.h" // Added  by Rushikesh

using namespace std;
using namespace astroaccelerate;

std::ofstream otf;

int main(int argc, const char *argv[])
{
  // new edits, part of working pipeline:(2)
  unsigned long int f_pos; // input file poistion variable

  int sysint, sysint1, sysint2, sysint3, sysint4, sysint5, sysint6; // this is to store value returned by system command when script is called at the end of file.

  // variables to read data from input file
  int device = 0;
  float sigma_cutoff, sigma_constant, max_boxcar_width_in_sec, periodicity_sigma_cutoff, periodicity_harmonics, max_dm_val;
  bool analysis_val = false, baselinenoise_val = false, debug_val = false;
  char file_name[128], config_path[128];

  // buf count loop variables
  int buf_count = 0;
  int new_file_flag = 1;

  // vaiable to enable storing of each output global_peaks.dat file
  int test_flg = 1;

  char cwd[256]; // to store current working directory

  std::string rm_cmnd, dir_cmnd, cat_cmnd;

  // timing variables
  float time_tmp, time_first_read, time_buf_read, time_buf_process, time_conc_files_find_frb, time_total;
  aa_gpu_timer timer_var, timer_var_total;

  float tsamp_unique;

  // variables to write the python file name
  std::stringstream stream;
  std::string s_ra, s_dec, s_mjd, s_bcount, py_fname, s_tsamp_unique, s_new_file_flag;

  // finding the intial precision to be set later
  std::streamsize ss = std::cout.precision();
  // std::cout << "Initial precision? = " << ss << '\n';

  aa_ddtr_plan ddtr_plan;

  // reading the input file for dm range and other input data::----------------------------------------------------
  // extracting the input file name from the argument passed through the terminal to the astroaccelerate.sh bash script
  if (argc < 1)
  {
    printf("\n\tUSAGE: ${ASTRO_ACCELERATE_EXECUTABLE_PATH}/examples_dedispersion_and_analysis ${input_file}\n");
    LOG(log_level::error, "Provide the AstroAccelerate configuration file. Exiting...\n");
    return -1;
  }

  // Reading the AstroAccelerate configuration file.
  aa_sigproc_input filterbank_data(argv[1]);
  filterbank_data.read_metadata_input(&device, &sigma_cutoff, &sigma_constant, &max_boxcar_width_in_sec, &periodicity_sigma_cutoff, &periodicity_harmonics, &ddtr_plan, &baselinenoise_val, &analysis_val, &debug_val, &file_name, &config_path);
  printf("\n ddtr_plan.user_dm(ddtr_plan.range()-1).high: %f \n", ddtr_plan.user_dm(ddtr_plan.range() - 1).high);

  aa_filterbank_metadata metadata = filterbank_data.read_metadata_from_FRB_SHM(); // FRB_SHM will be initialized for reading in this function.

  aa_device_info selected_device(device);

  //-------------- Configure pipeline. Select components and their options
  aa_pipeline::pipeline pipeline_components;
  pipeline_components.insert(aa_pipeline::component::dedispersion); // pipeline must always contain dedispersion step
  pipeline_components.insert(aa_pipeline::component::analysis);     // optional
  // pipeline_components.insert(aa_pipeline::component::periodicity); // optional
  // pipeline_components.insert(aa_pipeline::component::fdas); // optional

  aa_pipeline::pipeline_option pipeline_options;
  pipeline_options.insert(aa_pipeline::component_option::msd_baseline_noise);
  //----------------------------------------------------------------------<

  const bool enable_MSD_outlier_rejection = true;
  aa_analysis_plan::selectable_candidate_algorithm candidate_algorithm = aa_analysis_plan::selectable_candidate_algorithm::peak_find;
  candidate_algorithm = aa_analysis_plan::selectable_candidate_algorithm::peak_filtering;

  aa_pipeline_api<unsigned short> pipeline_runner(pipeline_components, pipeline_options, metadata, filterbank_data.input_buffer().data(), selected_device);

  pipeline_runner.bind(ddtr_plan);

  aa_analysis_plan analysis_plan(pipeline_runner.ddtr_strategy(), sigma_cutoff, sigma_constant, max_boxcar_width_in_sec, candidate_algorithm, enable_MSD_outlier_rejection);
  pipeline_runner.bind(analysis_plan);

  if (pipeline_runner.ready())
  {
    LOG(log_level::notice, "Pipeline is ready.");
  }
  else
  {
    LOG(log_level::notice, "Pipeline is not ready.");
    return 0; // new edits
  }

  // reading first buffer:

  timer_var_total.Start(); // starting timer for total time taken

  timer_var.Start(); // starting timer for first read

  //---- File name for raw file for 1K@3ms ------Rushikesh edits
  stream << std::fixed << std::setprecision(6) << metadata.src_raj();
  s_ra = stream.str();
  stream.str(std::string());
  stream.clear();

  stream << std::fixed << std::setprecision(6) << metadata.src_dej();
  s_dec = stream.str();
  stream.str(std::string());
  stream.clear();

  stream << std::fixed << std::setprecision(ss) << metadata.tstart();
  s_mjd = stream.str();
  stream.str(std::string());
  stream.clear();

  system("pwd > ./pwd.txt");
  ifstream ifs("./pwd.txt");
  std::string line, path;

  while (getline(ifs, line))
    path += line;
  printf("Output file path :: %s \n", path.c_str());

  std::string FileNameRaw;
  std::string FilePath;
  FileNameRaw = "Rawfile_1K@3ms_" + s_mjd + "_" + s_ra + "_" + s_dec + ".raw";
  FilePath = path + "/" + FileNameRaw;
  printf("##### FilePath ::: %s", FilePath.c_str());

  //------End------

  if (!filterbank_data.read_buf_count(buf_count, metadata))
  {
    std::cout << "ERROR: Could not read telescope data." << std::endl;
    return 0;
  }
  timer_var.Stop();

  time_first_read = timer_var.Elapsed() / 1000;
  printf("\n\n === First buffer read ===\n");
  //// RD: Check ////
  // NOT USED THIS ONE IN 2021 VERSION
  // aa_permitted_pipelines_2<aa_pipeline::component_option::zero_dm, false> runner(pipeline_runner.ddtr_strategy(), analysis_plan, filterbank_data.input_buffer().data());
  // IN 2021 VERSION analysis_plan REPALCED BY pipeline_runner.analysis_strategy()
  aa_permitted_pipelines_2<aa_pipeline::component_option::zero_dm, false> runner(pipeline_runner.ddtr_strategy(), pipeline_runner.analysis_strategy(), filterbank_data.input_buffer().data());

  //  bool dump_to_disk = false;
  bool dump_to_disk = true;
  bool dump_to_user = true;
  std::vector<analysis_output> output;

  //---------------------------------------------------------------------------------------------------------------------

  if (runner.setup())
  {

    for (buf_count = 0; true; buf_count++)
    {

      printf("\n buf count:: :: %d \n", buf_count);

      // read the timestamp of radec.dat
      if (buf_count == 0)
      {
        new_file_flag = 0;
      }

      // removing files from previous run and creating directories
      // if test_flg is true, it will create a bcount$ directory for each buffer ($ = buf_count)
      if (test_flg)
      {
        stream << std::setfill('0') << std::setw(4) << buf_count;
        s_bcount = stream.str();
        stream.str(std::string());
        stream.clear();

        rm_cmnd = "rm ./bcount" + s_bcount + "/global*";
        dir_cmnd = "bcount" + s_bcount;
        printf("\n buf count remove command: %s \n", rm_cmnd.c_str());
        printf("\n buf count mkdir command: %s \n", dir_cmnd.c_str());
        sysint2 = system("pwd");
        mkdir(dir_cmnd.c_str(), ACCESSPERMS);
        sysint2 = system(rm_cmnd.c_str());
      }

      sysint1 = system("rm -f analysed*");
      sysint2 = system("rm -f acc*");
      sysint3 = system("rm -f global*");
      sysint4 = system("rm -f fourier*");
      sysint5 = system("rm -f harmonic*");
      sysint6 = system("rm -f candidate*");
      sysint = system("rm -f peak*");
      printf("\n all files removed \n");

      // read only if the buffer number is greater than 0 - because first buffer has been read above already:
      if (buf_count > 0)
      {
        printf("\n we are here with buf_count: %d \n", buf_count);

        timer_var.Start();
        filterbank_data.open();
        if (!filterbank_data.read_buf_count(buf_count, metadata))
        {
          std::cout << "ERROR: Could not read telescope data." << std::endl;
          return 0;
        }

        // NOT USED THIS ONE IN 2021 VERSION
        // aa_permitted_pipelines_2<aa_pipeline::component_option::zero_dm, false> runner(pipeline_runner.ddtr_strategy(), analysis_plan, filterbank_data.input_buffer().data());
        //  IN 2021 VERSION analysis_plan REPALCED BY pipeline_runner.analysis_strategy()
        aa_permitted_pipelines_2<aa_pipeline::component_option::zero_dm, false> runner(pipeline_runner.ddtr_strategy(), pipeline_runner.analysis_strategy(), filterbank_data.input_buffer().data());

        timer_var.Stop();
        time_buf_read = timer_var.Elapsed() / 1000;
        printf("\n\n === Next buffer read ===\n");
      } // RD: 'buf_count > 0' loop ends here

      timer_var.Start();

      while (runner.next(dump_to_disk, dump_to_user, output))
      {
        LOG(log_level::notice, "Pipeline running over next chunk.");

        for (size_t i = 0; i < output.size(); i++)
        {
          std::cout << " dm_low " << output.at(i).dm_low << " dm_high " << output.at(i).dm_high << std::endl;
          for (size_t j = 0; j < output.at(i).pulses.size(); j++)
          {
            analysis_pulse p = output.at(i).pulses.at(j);
            std::cout << "dispersion measure " << p.dispersion_measure << " time " << p.time << " snr " << p.snr << " pulse_width " << p.pulse_width << std::endl;
          }
        }
      }

      timer_var.Stop();
      time_buf_process = timer_var.Elapsed() / 1000;
      printf("\n\n === Buffer processed ===\n");

      // creating python output file name
      printf("\n RA: %f DEC: %f MJD: %f \n", metadata.src_raj(), metadata.src_dej(), metadata.tstart());

      stream << new_file_flag;
      s_new_file_flag = stream.str();
      stream.str(std::string());
      stream.clear();

      stream << std::setfill('0') << std::setw(4) << buf_count;
      s_bcount = stream.str();
      stream.str(std::string());
      stream.clear();

      stream << std::fixed << std::setprecision(6) << metadata.src_raj();
      s_ra = stream.str();
      stream.str(std::string());
      stream.clear();

      stream << std::fixed << std::setprecision(6) << metadata.src_dej();
      s_dec = stream.str();
      stream.str(std::string());
      stream.clear();

      stream << std::fixed << std::setprecision(ss) << metadata.tstart();
      s_mjd = stream.str();
      stream.str(std::string());
      stream.clear();

      stream << std::fixed << std::setprecision(ss) << tsamp_unique;
      s_tsamp_unique = stream.str();
      stream.str(std::string());
      stream.clear();

      printf("\n After setting precision and string conversion: RA: %s DEC: %s MJD: %s buf_count: %s \n", s_ra.c_str(), s_dec.c_str(), s_mjd.c_str(), s_bcount.c_str());

      py_fname = "python "; // --bc 0 --ra 2.5 --dec -36 --mjd 55.5";
      // py_fname += configpath_name;
      py_fname += "spsplotii_frb_detect.py global_peaks.dat";
      py_fname += " --mjd " + s_mjd + " --ra " + s_ra + " --dec " + s_dec + " --bc " + s_bcount + " --buf_t " + s_tsamp_unique + " --new_file_flag " + s_new_file_flag;
      printf("\n py_fname: %s \n", py_fname.c_str());

      // writing output files and running python script
      timer_var.Start();

      sysint1 = system("cat peak* > global_peaks.dat");
      sysint2 = system("cat analysed* > global_analysed_frb.dat");
      sysint3 = system("cat fourier-* > global_periods.dat");
      sysint4 = system("cat fourier_inter* > global_interbin.dat");
      sysint5 = system("cat harmo* > global_harmonics.dat");
      sysint6 = system("cat candidate* > global_candidates.dat");
      sysint = system(py_fname.c_str());

      /*-------RD: To write python output file -------
           const char * ls_args[2] = { py_fname.c_str(), NULL} ;
           pid_t c_pid, pid;
           int status;

           c_pid = fork();

           if (c_pid == 0){
             /* CHILD */

      /*       printf("Child: executing ls\n");
             execvp( ls_args[0], ls_args);

             perror("execve failed");
             }else if (c_pid > 0){
              /* PARENT */

      /*        if( (pid = wait(&status)) < 0){
                perror("wait");
                _exit(1);
              }

              printf("Parent: finished\n");

            }else{
              perror("fork failed");
              _exit(1);
            }

      //--------End------------------------*/

      printf("\nOutput files written\n");

      // if test_flg is true, it will save every global_peaks.dat file from every buffer
      if (test_flg)
      {
        cat_cmnd = "cat peak* > ./bcount" + s_bcount + "/global_peaks.dat";
        printf("\n Concatenation command:: %s \n", cat_cmnd.c_str());
        sysint2 = system(cat_cmnd.c_str());
        int alp, prod, block_no;
        alp = stoi(s_bcount);
        prod = alp * 7;
        if (prod > 16)
          block_no = prod % 16;
        else
          block_no = prod;
        // memcpy(&timestamp_gps, &(dataHdr->timestamp_gps[block_no]), sizeof(timestamp_gps));
        // string filepath;
        // time_t time_l = timestamp_gps.tv_sec;
        // char *dt = ctime(&time_l);
        // filepath = "./bcount" + s_bcount + "/TimeStamp.txt";
        // std::ofstream outputFile(filepath.c_str());
        // outputFile << "Timestamp:\t" << dt << timestamp_gps.tv_usec << " microseconds.\n";
        // outputFile.close();
      }

      timer_var.Stop();
      time_conc_files_find_frb = timer_var.Elapsed() / 1000;
      printf("\n\n === Output files created | Concatenated | Python script ran ===\n");

      if (buf_count == 0)
      {
        otf.open("/home/frbuser/Viren/testf2.txt", std::ios_base::app);
        otf << time_first_read << endl;
        otf.close();
        otf.open("/home/frbuser/Viren/testp2.txt", std::ios_base::app);
        otf << time_buf_process << endl;
        otf.close();
        printf("\n--------------------- Timing information after processing first buffer ---------------------\n");
        printf("\n Time to read first buffer: %f (GPU estimate) \n", time_first_read);
        printf("\n Time to process buffer number %d: %f (GPU estimate) \n", buf_count, time_buf_process);
        printf("\n Time to concatenate files and run python script for buffer number %d: %f (GPU estimate) \n", buf_count, time_conc_files_find_frb);
        printf("\n--------------------------------------------------------------------------------------------\n");
      }
      else
      {
        otf.open("/home/frbuser/Viren/testf2.txt", std::ios_base::app);
        otf << time_buf_read << endl;
        otf.close();
        otf.open("/home/frbuser/Viren/testp2.txt", std::ios_base::app);
        otf << time_buf_process << endl;
        otf.close();
        printf("\n--------------------- Timing information after processing new buffer ---------------------\n");
        printf("\n Time to read buffer number %d: %f (GPU estimate) \n", buf_count, time_buf_read);
        printf("\n Time to process buffer number %d: %f (GPU estimate) \n", buf_count, time_buf_process);
        printf("\n Time to concatenate files and run python script for buffer number %d: %f (GPU estimate) \n", buf_count, time_conc_files_find_frb);
        printf("\n-------------------------------------------------------------------------------------------\n");
      }
    }
  }

  timer_var_total.Stop();
  time_total = timer_var_total.Elapsed() / 1000;

  printf("\n------------------------------- Final Timing Information  ---------------------------------\n");
  printf("\n Total number of buffers processed: %d \n", buf_count);
  printf("\n Total time taken: %f (GPU estimate) \n", time_total);
  printf("\n-------------------------------------------------------------------------------------------\n");

  LOG(log_level::notice, "Finished.");
  return 0;
}