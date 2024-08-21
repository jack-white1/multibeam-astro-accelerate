#include "aa_sigproc_input.hpp"
#include "aa_log.hpp"
#include "aa_ddtr_plan.hpp" // added by abhinav
#include "externalLibraries.h"
#include "frb_shm.h"
#include "gmrt_newcorr.h"

namespace astroaccelerate
{
  BeamHeaderType *dataHdr_FRB;
  DataBuffer *dataBuffer_FRB;
  int RecNum_FRB = 0;
  unsigned int DataSeq_FRB = 0;
  std::vector<unsigned char> rawData(DataSize);

  /** \brief Constructor for aa_sigproc_input. */
  aa_sigproc_input::aa_sigproc_input(const std::string &path) : header_is_read(false), data_is_read(false), file_path(path)
  {
    isopen = false;
  }

  /** \brief Destructor for aa_sigproc_input. */
  aa_sigproc_input::~aa_sigproc_input()
  {
    if (isopen)
    {
      close();
    }
  }

  int aa_sigproc_input::initialize_FRB_SHM()
  {
    int idDataHdr, idDataBuffer;
    idDataHdr = shmget(DasHeaderKey, sizeof(BeamHeaderType), SHM_RDONLY);
    if (idDataHdr < 0)
    {
      perror("FRB shmget HEADER ID");
      return -1;
    }
    dataHdr_FRB = (BeamHeaderType *)shmat(idDataHdr, 0, SHM_RDONLY);
    if ((void *)dataBuffer_FRB == (void *)-1)
    {
      perror("FRB shmat HEADER ID");
      return -3;
    }

    idDataBuffer = shmget(DasBufferKey, sizeof(DataBuffer), SHM_RDONLY);
    if (idDataHdr < 0)
    {
      perror("FRB shmget BUFFER ID");
      return -1;
    }
    dataBuffer_FRB = (DataBuffer *)shmat(idDataBuffer, 0, SHM_RDONLY);
    if ((void *)dataBuffer_FRB == (void *)-1)
    {
      perror("FRB shmat BUFFER ID");
      return -3;
    }

    if (dataBuffer_FRB->is_buf_empty)
    {
      RecNum_FRB = dataBuffer_FRB->curRecord;
      DataSeq_FRB = dataBuffer_FRB->curBlock;
    }
    else
    {
      RecNum_FRB = (dataBuffer_FRB->curRecord - 1) % MaxDataBlocks;
      DataSeq_FRB = dataBuffer_FRB->curBlock - 1;
    }

    m_input_buffer.resize((size_t)DataSize * (size_t)Blocks_in_Buf_count);

    LOG(log_level::notice, "RecNum_FRB = " + std::to_string(RecNum_FRB) + "; dataBuffer_FRB->curRecord = " + std::to_string(dataBuffer_FRB->curRecord));
    LOG(log_level::notice, "DataSeq_FRB = " + std::to_string(DataSeq_FRB) + "; dataBuffer_FRB->curBlock = " + std::to_string(dataBuffer_FRB->curBlock));

    return 0;
  }

  /**
   * \brief Method to open the sigproc input file.
   * \returns A boolean flag to indicate whether the operation was successful (true) or not (false).
   */
  bool aa_sigproc_input::open()
  {
    fp = fopen(file_path.c_str(), "rb");
    if (fp == NULL)
    {
      return false;
    }
    isopen = true;
    return true;
  }

  /**
   * \brief Closes the input file.
   * \returns A boolean flag to indicate whether the operation was successful (true) or not (false).
   */
  bool aa_sigproc_input::close()
  {
    if (isopen)
    {
      if (fclose(fp) == -1)
      {
        return false;
      }
      isopen = false;
    }
    return true;
  }

  /** \brief Method to read the AstroAccelerate configuration from GMRT.SPS.10m.newAA.txt file. **/
  bool aa_sigproc_input::read_metadata_input(int *selected_card_id, float *sigma_cutoff_input, float *sigma_constant_input, float *max_boxcar_width_in_sec_input, float *periodicity_sigma_cutoff_input, float *periodicity_harmonics_input, aa_ddtr_plan *ddtr_plan, bool *baselinenoise_val, bool *analysis_val, bool *debug_val, char (*file_name_input)[128], char (*config_path)[128])
  {
    if (!isopen)
    {
      if (!open())
      {
        printf("\n returning empty handed \n");
        return false;
      }
    }

    char line[128];
    aa_ddtr_plan d_p_input;
    float high_input, low_input, step_input;
    int inBin_input, outBin_input;

    while (fgets(line, sizeof(line), fp))
    {
      if (strncmp(line, "range", 5) == 0) // Frequency of channel 1 (MHz)
      {
        low_input = strtof(&line[9], NULL);
        high_input = strtof(&line[14], NULL);
        step_input = strtof(&line[19], NULL);
        inBin_input = strtol(&line[24], NULL, 10);
        outBin_input = strtol(&line[26], NULL, 10);
        printf("\n low : %f high : %f step : %f inBin : %d outBin : %d \n", low_input, high_input, step_input, inBin_input, outBin_input);
        ddtr_plan->add_dm(low_input, high_input, step_input, inBin_input, outBin_input); // this use of -> is quite a find !! abhinav
      }
      else if (strncmp(line, "sigma_cutoff", 12) == 0) // Frequency of channel 1 (MHz)
      {
        *sigma_cutoff_input = strtof(&line[13], NULL); // new edits, was (&line[16], NULL) before
        printf("\n sigma_cutoff :  ::: %f \n", *sigma_cutoff_input);
      }
      else if (strncmp(line, "sigma_constant", 14) == 0) // Time stamp of first sample (MJD)
      {
        *sigma_constant_input = strtof(&line[16], NULL);
        printf("\n sigma_constant_input :  ::: %f \n", *sigma_constant_input);
      }
      else if (strncmp(line, "max_boxcar_width_in_sec", 23) == 0) // "Channel bandwidth"
      {
        *max_boxcar_width_in_sec_input = strtof(&line[24], NULL);
        printf("\n max_boxcar_width_in_sec : :::%f \n", *max_boxcar_width_in_sec_input);
      }
      else if (strncmp(line, "periodicity_sigma_cutoff", 24) == 0) // "Channel bandwidth"
      {
        *periodicity_sigma_cutoff_input = strtof(&line[25], NULL);
        printf("\n periodicity_sigma_cutoff_input: ::: %f \n", *periodicity_sigma_cutoff_input);
      }
      else if (strncmp(line, "periodicity_harmonics", 21) == 0) // "Channel bandwidth"
      {
        *periodicity_harmonics_input = strtof(&line[22], NULL);
        printf("\n periodicity_harmonics_input: ::: %f \n", *periodicity_harmonics_input);
      }
      else if (strncmp(line, "baselinenoise", 13) == 0) // "Channel bandwidth"
      {
        *baselinenoise_val = true;
        printf("\n baselinenoise_val: ::: %d \n", *baselinenoise_val);
      }
      else if (strncmp(line, "file", 4) == 0) // "Channel bandwidth"
      {
        strcpy(*file_name_input, &line[5]);
        printf("\n line is: ::: %s :: file_name_input is :: %s :: \n", line, *file_name_input);
      }
      else if (strncmp(line, "config_file_path", 16) == 0) // "Channel bandwidth"
      {
        strcpy(*config_path, &line[17]);
        // printf("\n config path line is: ::: %s \n", line);
      }
    }
    return true;
  }

  aa_filterbank_metadata aa_sigproc_input::read_metadata_from_FRB_SHM()
  {
    LOG(log_level::notice, "Initializing the FRB_SHM.");

    initialize_FRB_SHM();

    LOG(log_level::notice, "Reading metadata from the FRB_SHM Header.");
    LOG(log_level::notice, "Update the variable ScanTab[0] with the correct index (instead of 0) based on the porject code. By SSK.\n");

    double az_start = 0;          // Useless for AstroAccelerate
    double za_start = 0;          // Useless for AstroAccelerate
    int BeamHostID;               // Unused as of July 12, 2024
    char BeamHostName[128];       // Unused as of July 12, 2024
    std::string source_name = ""; // Unused as of July 12, 2024
    unsigned int antmask = 0;     // Unused as of July 12, 2024
    double src_raj = 0;           // Unused as of July 12, 2024
    double src_dej = 0;           // Unused as of July 12, 2024
    double tstart = 0;            // timestamp (MJD) of the first sample
    double tsamp = 0.00131072;    // time interval between samples
    double refdm = 0;             // Useless for AstroAccelerate
    double period = 0;            // Useless for AstroAccelerate
    double fch1 = 0;              // in MHz
    double foff = 0;              // in MHz
    float bandwidth = 0;          // in MHz // Unused as of July 12, 2024
    double sign_bw = 0;           // RD : inv_check // From Hdr SHM
    double fchannel = 0;          // Useless for AstroAccelerate
    int telescope_id = 7;         // GMRT has telescope_id = 7. Refer https://github.com/lbaehren/lofarsoft/blob/master/src/Pulsar/external/dspsr/dspsr-20120327/Kernel/Formats/sigproc/SigProcObservation.C
    int machine_id = 0;           // 0=FAKE
    int data_type = 1;            // Useless for AstroAccelerate
    int barycentric = 0;          // Useless for AstroAccelerate
    int pulsarcentric = 0;        // Useless for AstroAccelerate
    int nbits = 8;
    int nsamples = ((int)FFT_Samps_per_Block * TEL_to_FRB_Block_factor * Blocks_in_Buf_count);
    int nchans = NCHANNELS;
    int nifs = 1;
    int strat = 0;   // Useless for AstroAccelerate
    float ovrlp = 1; // Useless for AstroAccelerate

    BeamHostID = dataHdr_FRB->BeamGenHdr.BeamHostID;
    *BeamHostName = *(dataHdr_FRB->BeamGenHdr.BeamHostName);
    // source_name = dataHdr_FRB->ScanTab[0].source.object;
    source_name = "B0329+54";
    // src_raj = dataHdr_FRB->ScanTab[0].source.ra_app;
    // src_dej = dataHdr_FRB->ScanTab[0].source.dec_app;
    // src_raj = dataHdr_FRB->BeamGenHdr.BeamSteeringParams.RA[0];  // Beam RA.
    // src_dej = dataHdr_FRB->BeamGenHdr.BeamSteeringParams.DEC[0]; // Beam Dec.
    src_raj = 53.24754;
    src_dej = 54.5787025;
    antmask = dataHdr_FRB->ScanTab[0].source.antmask;
    tstart = 60513.93; // 22 July 2024, 22:18 hrs
    // tstart = dataBuffer_FRB->timestamp_gps->tv_sec + (dataBuffer_FRB->timestamp_gps->tv_usec / 1e6) + (dataBuffer_FRB->blk_nano[dataBuffer_FRB->curRecord] / 1e9);
    // tsamp = dataHdr_FRB->corr.daspar.gsb_final_bw * dataHdr_FRB->BeamGenHdr.SampInterval * dataHdr_FRB->BeamGenHdr.PostTimeInt[BeamHostID] / (dataHdr_FRB->corr.corrpar.clock / 1e6);
    // bandwidth = dataHdr_FRB->scanpar->BW;
    // bandwidth = dataHdr_FRB->corr.daspar.gsb_acq_bw / dataHdr_FRB->corr.daspar.gsb_final_bw;
    /* Comment the next three lines*/
    fch1 = 500;
    sign_bw = 1;
    foff = 200 / NCHANNELS;
    /* Uncomment the next three lines*/
    // fch1 = dataHdr_FRB->ScanTab[0].source.freq[0];
    // sign_bw = -1 * (dataHdr_FRB->ScanTab[0].source.net_sign[0]);
    // foff = dataHdr_FRB->ScanTab[0].source.ch_width * sign_bw;

    aa_filterbank_metadata meta(telescope_id,
                                machine_id,
                                data_type,
                                "unknown", /*rawdatafile */
                                source_name,
                                barycentric,
                                pulsarcentric,
                                az_start,
                                za_start,
                                src_raj,
                                src_dej,
                                tstart,
                                tsamp,
                                nbits,
                                nsamples,
                                fch1,
                                foff,
                                sign_bw,
                                0, /*FREQUENCY_START */
                                fchannel,
                                0, /*FREQUENCY_END */
                                nchans,
                                nifs,
                                refdm,
                                period,
                                strat,
                                ovrlp);

    header_is_read = true;
    return meta;
  }

  /**
   * \brief If the file is open, and the header has been read, and the data has not yet been read, then read the input data from the input file.
   * \details Reading the telescope input data can only be performed once, after which this method will always return false.
   * \returns A boolean flag to indicate whether the operation was successful (true) or not (false).
   * \warning The method will return true only once, that is the first time the data are read from the input successfully. At this point the input_buffer should be checked for data.
   */
  bool aa_sigproc_input::read_signal()
  {
    /*    if(!isopen || !header_is_read || data_is_read) {
          return false;
        }
    */
    // new edit(check)

    get_recorded_data(m_input_buffer);
    return true;
  }

  template <typename T>
  bool aa_sigproc_input::get_recorded_data(std::vector<T> &input_buffer)
  {
    const size_t inputsize = (size_t)m_meta.nsamples() * (size_t)m_meta.nchans();
    input_buffer.resize(inputsize);
    int c;

    unsigned long int total_data;
    const int nchans = m_meta.nchans();
    const int nbits = m_meta.nbits();
    //{{{ Load in the raw data from the input file and transpose
    if (nbits == 32)
    {
      // Allocate a tempory buffer to store a line of frequency data
      float *temp_buffer = (float *)malloc(nchans * sizeof(float));

      // Allocate floats to hold the mimimum and maximum value in the input data
      float max = -10000000000;
      float min = 10000000000;

      // Allocate a variable to hold the file pointer position
      fpos_t file_loc;

      // Set the file pointer position
      fgetpos(fp, &file_loc);

      // Find the minimum and maximum values in the input file.
      while (!feof(fp))
      {
        if (fread(temp_buffer, sizeof(float), nchans, fp) != (size_t)nchans)
          break;
        for (c = 0; c < nchans; c++)
        {
          if (temp_buffer[c] > max)
            max = temp_buffer[c];
          if (temp_buffer[c] < min)
            min = temp_buffer[c];
        }
      }

      // Calculate the bin size in a distribution of unsigned shorts for the input data.
      float bin = (max - min) / 65535.0f;

      printf("\n Conversion Parameters: %f\t%f\t%f", min, max, bin);

      // Move the file pointer back to the end of the header
      fsetpos(fp, &file_loc);

      // Read in the data, scale to an unsigned short range, transpose it and store it in the input buffer
      total_data = 0;
      while (!feof(fp))
      {
        if (total_data % (int)(m_meta.nsamples() * 0.1) == 0)
          printf("Reading record %10lu / %10d \t (%2d%% complete)\n", total_data, m_meta.nsamples(), (int)(100.0 * ((float)total_data / (float)m_meta.nsamples())));
        if (fread(temp_buffer, sizeof(float), nchans, fp) != (size_t)nchans)
          break;
        for (c = 0; c < nchans; c++)
        {
          (input_buffer)[c + total_data * (nchans)] = (unsigned short)((temp_buffer[c] - min) / bin);
        }
        total_data++;
      }
      printf("Reading record %10lu / %10d \t (%2d%% complete)\n", total_data, m_meta.nsamples(), (int)(100.0 * ((float)total_data / (float)m_meta.nsamples())));
      free(temp_buffer);
    }
    else if (nbits == 16)
    {
      // Allocate a tempory buffer to store a line of frequency data
      unsigned short *temp_buffer = (unsigned short *)malloc(nchans * sizeof(unsigned short));

      // Read in the data, transpose it and store it in the input buffer
      total_data = 0;
      while (!feof(fp))
      {
        if (total_data % (int)(m_meta.nsamples() * 0.1) == 0)
          printf("Reading record %10lu / %10d \t (%2d%% complete)\n", total_data, m_meta.nsamples(), (int)(100.0 * ((float)total_data / (float)m_meta.nsamples())));
        if (fread(temp_buffer, sizeof(unsigned short), nchans, fp) != (size_t)nchans)
          break;
        for (c = 0; c < nchans; c++)
        {
          (input_buffer)[c + total_data * (nchans)] = (unsigned short)temp_buffer[c];
        }
        total_data++;
      }
      printf("Reading record %10lu / %10d \t (%2d%% complete)\n", total_data, m_meta.nsamples(), (int)(100.0 * ((float)total_data / (float)m_meta.nsamples())));
      free(temp_buffer);
    }
    else if (nbits == 8)
    {
      /* new edits(original code)
            // Allocate a tempory buffer to store a line of frequency data
            unsigned char *temp_buffer = (unsigned char *) malloc(nchans * sizeof(unsigned char));

            // Read in the data, transpose it and store it in the input buffer
            total_data = 0;
            while (!feof(fp)) {
        if(total_data % (int)(m_meta.nsamples()*0.1) == 0) printf("Reading record %10lu / %10d \t (%2d%% complete)\n", total_data, m_meta.nsamples(), (int)(100.0*((float)total_data/(float)m_meta.nsamples())));
        if (fread(temp_buffer, sizeof(unsigned char), nchans, fp) != (size_t)nchans) {
          break;
        }
        for (c = 0; c < nchans; c++) {
          ( input_buffer )[c + total_data * ( nchans )] = (unsigned short) temp_buffer[c];
        }
        total_data++;
            }
            printf("Reading record %10lu / %10d \t (%2d%% complete)\n", total_data, m_meta.nsamples(), (int)(100.0*((float)total_data/(float)m_meta.nsamples())));
            free(temp_buffer);
          }
      */

      // new edits

      printf("\n We are here with nbits: %d\n", nbits);

      fseek(fp, 399, SEEK_SET);

      printf("\n Setting pointer to 399.\n");

      // Allocate a tempory buffer to store a line of frequency data
      unsigned char *temp_buffer = (unsigned char *)malloc(nchans * sizeof(unsigned char));

      // Read in the data, transpose it and store it in the input buffer
      total_data = 0;

      while (!feof(fp))
      {

        //      while (!feof(fp)) {
        if (total_data % (int)(m_meta.nsamples() * 0.1) == 0)
          printf("Reading record %10lu / %10d \t (%2d%% complete)\n", total_data, m_meta.nsamples(), (int)(100.0 * ((float)total_data / (float)m_meta.nsamples())));
        //	if(total_data % (int)(nsamples*0.1) == 0) printf("Reading record %10lu / %10d \t (%2d%% complete)\n", total_data, nsamples, (int)(100.0*((float)total_data/(float)nsamples)));

        if (fread(temp_buffer, sizeof(unsigned char), nchans, fp) != (size_t)nchans)
        {

          break;
        }

        for (c = 0; c < nchans; c++)
        {

          (input_buffer)[c + total_data * (nchans)] = (unsigned short)temp_buffer[c];
          //          if(c < 700 && total_data == 0) printf("\n %hu \n", ( input_buffer )[c + total_data * ( nchans )]);
        }
        total_data++;

/*        if (total_data == nsamples) 
         {
           break;
         }
*/      }
free(temp_buffer);
    }

    // end
    else if (nbits == 4)
    {
      // Allocate a temporary buffer to store a line of frequency data
      // each byte stores 2 frequency data
      int nb_bytes = nchans / 2;
      unsigned char *temp_buffer = (unsigned char *)malloc(nb_bytes * sizeof(unsigned char));
      // Read in the data, transpose it and store it in the input buffer
      total_data = 0;
      // 00001111
      char mask = 0x0f;
      while (!feof(fp))
      {
        if (total_data % (int)(m_meta.nsamples() * 0.1) == 0)
          printf("Reading record %10lu / %10d \t (%2d%% complete)\n", total_data, m_meta.nsamples(), (int)(100.0 * ((float)total_data / (float)m_meta.nsamples())));
        if (fread(temp_buffer, sizeof(unsigned char), nb_bytes, fp) != (size_t)nb_bytes)
          break;
        for (c = 0; c < nb_bytes; c++)
        {
          // (n >> a) & ( (1 << a) - 1) -> right shift by 'a' bits, then keep the last 'b' bits
          // Here, we right shift 4 bits and keep the last 4 bits
          (input_buffer)[(c * 2) + total_data * (nchans)] = (unsigned short)((temp_buffer[c] >> 4) & mask);
          // n & ( (1 << a ) - 1)
          // Here we keep the last 4 bits
          (input_buffer)[(c * 2) + 1 + total_data * (nchans)] = (unsigned short)(temp_buffer[c] & mask);
        }
        total_data++;
      }
      printf("Reading record %10lu / %10d \t (%2d%% complete)\n", total_data, m_meta.nsamples(), (int)(100.0 * ((float)total_data / (float)m_meta.nsamples())));
      free(temp_buffer);
    }
    else if (nbits == 2)
    {
      // Allocate a temporary buffer to store a line of frequency data
      // each byte stores 4 frequency data
      int nb_bytes = nchans / 4;
      unsigned char *temp_buffer = (unsigned char *)malloc(nb_bytes * sizeof(unsigned char));
      // Read in the data, transpose it and store it in the input buffer
      total_data = 0;
      // 00001111
      char mask = 0x03;
      while (!feof(fp))
      {
        if (total_data % (int)(m_meta.nsamples() * 0.1) == 0)
          printf("Reading record %10lu / %10d \t (%2d%% complete)\n", total_data, m_meta.nsamples(), (int)(100.0 * ((float)total_data / (float)m_meta.nsamples())));
        if (fread(temp_buffer, sizeof(unsigned char), nb_bytes, fp) != (size_t)nb_bytes)
          break;
        for (c = 0; c < nb_bytes; c++)
        {
          (input_buffer)[(c * 4) + total_data * (nchans)] = (unsigned short)((temp_buffer[c] >> 6) & mask);
          (input_buffer)[(c * 4) + 1 + total_data * (nchans)] = (unsigned short)((temp_buffer[c] >> 4) & mask);
          (input_buffer)[(c * 4) + 2 + total_data * (nchans)] = (unsigned short)((temp_buffer[c] >> 2) & mask);
          (input_buffer)[(c * 4) + 3 + total_data * (nchans)] = (unsigned short)(temp_buffer[c] & mask);
        }
        total_data++;
      }
      printf("Reading record %10lu / %10d \t (%2d%% complete)\n", total_data, m_meta.nsamples(), (int)(100.0 * ((float)total_data / (float)m_meta.nsamples())));
      free(temp_buffer);
    }
    else if (nbits == 1)
    {
      // each byte stores 8 frequency data
      int nb_bytes = nchans / 8;

      // Allocate a temporary buffer to store a line of frequency data
      unsigned char *temp_buffer = (unsigned char *)malloc(nb_bytes * sizeof(unsigned char));

      // Read in the data, transpose it and store it in the input buffer
      total_data = 0;

      while (!feof(fp))
      {
        if (total_data % (int)(m_meta.nsamples() * 0.1) == 0)
          printf("Reading record %10lu / %10d \t (%2d%% complete)\n", total_data, m_meta.nsamples(), (int)(100.0 * ((float)total_data / (float)m_meta.nsamples())));
        if (fread(temp_buffer, sizeof(unsigned char), nb_bytes, fp) != (size_t)nb_bytes)
          break;

        for (c = 0; c < nb_bytes; c++)
        {
          for (int i = 0; i < 8; i++)
          {
            unsigned char mask = 1 << i;
            unsigned char masked_char = temp_buffer[c] & mask;
            unsigned char bit = masked_char >> i;
            (input_buffer)[(c * 8) + i + total_data * (nchans)] = (unsigned short)bit;
          }
        }
        total_data++;
      }
      printf("Reading record %10lu / %10d \t (%2d%% complete)\n", total_data, m_meta.nsamples(), (int)(100.0 * ((float)total_data / (float)m_meta.nsamples())));
      free(temp_buffer);
    }
    else
    {
      printf("ERROR: Invalid number of bits in input data.\n");
      return false;
    }

    data_is_read = true;
    return true;
  }

  bool aa_sigproc_input::read_buf_count(int buf_count, const aa_filterbank_metadata &filterbank_metadata)
  {
    LOG(log_level::notice, "\n------------------ Reading buffer count: " + std::to_string(buf_count) + " ------------------\n");

    unsigned char *buffer = (unsigned char *)malloc(NCHANNELS);
    unsigned char temp, *blkp_FRB;
    int iSamp, iCh, beam_ID = 0, overlap = 0;
    long i = 0;
    if (buf_count != 0)
      overlap = 1;

    while (dataBuffer_FRB->is_buf_empty)
    {
      sleep(2);
      LOG(log_level::notice, "Sleeping peacefully since the buffer is full of garbage.");
      RecNum_FRB = dataBuffer_FRB->curRecord;
      DataSeq_FRB = dataBuffer_FRB->curBlock;
    }

    RecNum_FRB = (dataBuffer_FRB->curRecord - 1 - overlap) % MaxDataBlocks;
    DataSeq_FRB = dataBuffer_FRB->curBlock - 1;

    for (int RecNum_in_bcnt = 0; RecNum_in_bcnt < Blocks_in_Buf_count; RecNum_in_bcnt++) // Read one buffer count in a go.
    {
      LOG(log_level::notice, "Buffer block being read: " + std::to_string(RecNum_in_bcnt));
      LOG(log_level::notice, "FRB_SHM block being read: " + std::to_string(RecNum_FRB));

      while (dataBuffer_FRB->active == 0)
      {
        LOG(log_level::notice, "FRB_SHM is not active, probably because the telescope is not observing any source of interest. Sleeping peacefully...");
        sleep(2); // Wait for 2 sec.
      }

      while (DataSeq_FRB >= dataBuffer_FRB->curBlock) // Waiting for the data to be written in the memory.
      {
        LOG(log_level::notice, "Sleeping peacefully...");
        sleep(2); // Wait for 2 sec.
      }
      if (dataBuffer_FRB->curBlock - DataSeq_FRB >= MaxDataBlocks - 1) // Realigning.
      {
        LOG(log_level::warning, "Processing lagged behind... Some buffer lost. Realigning to the latest buffer position.");
        RecNum_FRB = (dataBuffer_FRB->curRecord - 1 + MaxDataBlocks) % MaxDataBlocks;
        DataSeq_FRB = dataBuffer_FRB->curBlock - 1;
      }

      blkp_FRB = dataBuffer_FRB->data + ((long)DataSize * NBeams * RecNum_FRB) + ((long)DataSize * beam_ID);
      memcpy(&rawData[0], blkp_FRB, (long)DataSize); // Refer to the FRB_SHM structure document.

      if (filterbank_metadata.sign_bw() == -1) // inverting the channels of the band if it is upper side-band
      {
        for (iSamp = 1; iSamp <= (int)FFT_Samps_per_Block; iSamp++)
        {
          memcpy(buffer, &rawData[0] + (iSamp - 1) * NCHANNELS, NCHANNELS);
          for (iCh = 0; iCh < NCHANNELS / 2; iCh++)
          {
            temp = buffer[iCh];
            buffer[iCh] = buffer[NCHANNELS - iCh - 1];
            buffer[NCHANNELS - iCh - 1] = temp;
          }
          memcpy(&rawData[0] + (iSamp - 1) * NCHANNELS, buffer, NCHANNELS);
        }
      }

      memcpy(&m_input_buffer[RecNum_in_bcnt * (long)DataSize], &rawData[0], DataSize);

      RecNum_FRB = (RecNum_FRB + 1) % MaxDataBlocks;
      DataSeq_FRB++;
    }

    free(buffer);
    data_is_read = true;
    return true;
  }
}