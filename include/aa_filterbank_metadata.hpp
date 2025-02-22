#ifndef ASTRO_ACCELERATE_AA_FILTERBANK_METADATA_HPP
#define ASTRO_ACCELERATE_AA_FILTERBANK_METADATA_HPP

#include <string>

#include "aa_log.hpp"

namespace astroaccelerate {

  /**
   * \class aa_filterbank_metadata aa_filterbank_metadata.hpp "include/aa_filterbank_metadata.hpp"
   * \brief Contains SIGPROC metadata.
   * \brief For data descriptions, please see SIGPROC interface control document.
   * \details Filterbank metadata must be provided for any permitted pipeline to run any component.
   * \details Data description from "SIGPROC-v3.7 (Pulsar) Signal Processing Programs".
   * \details Source: http://sigproc.sourceforge.net/sigproc.pdf.
   * \author Cees Carels.
   * \date 5 November 2018.
   */

  class aa_filterbank_metadata {
  public:

    /** \brief Trivial constructor for aa_filterbank_metadata. */
    aa_filterbank_metadata() : m_az_start(0),
			       m_za_start(0),
			       m_src_raj(0),
			       m_src_dej(0),
			       m_tstart(0),
			       m_tsamp(0),
			       m_refdm(0),
			       m_period(0),
			       m_fch1(0),
			       m_foff(0),
			       m_sign_bw(0), // RD: inv_check //new edits
			       m_fchannel(0),
			       m_telescope_id(0),
			       m_machine_id(0),
			       m_data_type(0),
			       m_barycentric(0),
			       m_pulsarcentric(0),
			       m_nbits(0),
			       m_nsamples(0),
			       m_nchans(0),
			       m_nifs(0),
			       m_FREQUENCY_START(0),
			       m_FREQUENCY_END(0),
			       m_rawdatafile(""),
			       m_source_name(""), 
			       m_strat(0), //new edits
			       m_ovrlp(0), //new edits
			       m_max_block_read(0),			//edits by Samyak
			       m_max_block_count(0),
			       m_tstamp_per_block(0),
			       m_block_ovrlp_band3(0),
			       m_block_ovrlp_band4(0),
			       m_block_ovrlp_band5(0),
			       m_enable_flag(0),
			       m_free_memory(0){  			// edits by Samyak end      
   }
  
    /** \brief Copy constructor for aa_filterbank_metadata. */
    aa_filterbank_metadata(const aa_filterbank_metadata &meta) : m_az_start(meta.m_az_start),
								 m_za_start(meta.m_za_start),
								 m_src_raj(meta.m_src_raj),
								 m_src_dej(meta.m_src_dej),
								 m_tstart(meta.m_tstart),
								 m_tsamp(meta.m_tsamp),
								 m_refdm(meta.m_refdm),
								 m_period(meta.m_period),
								 m_fch1(meta.m_fch1),
								 m_foff(meta.m_foff),
								 m_sign_bw(meta.m_sign_bw), //new edits
								 m_fchannel(meta.m_fchannel),
								 m_telescope_id(meta.m_telescope_id),
								 m_machine_id(meta.m_machine_id),
								 m_data_type(meta.m_data_type),
								 m_barycentric(meta.m_barycentric),
								 m_pulsarcentric(meta.m_pulsarcentric),
								 m_nbits(meta.m_nbits),
								 m_nsamples(meta.m_nsamples),
								 m_nchans(meta.m_nchans),
								 m_nifs(meta.m_nifs),
								 m_FREQUENCY_START(meta.m_FREQUENCY_START),
								 m_FREQUENCY_END(meta.m_FREQUENCY_END),
								 m_rawdatafile(meta.m_rawdatafile),
								 m_source_name(meta.m_source_name),
								 m_strat(meta.m_strat), //new edits
								 m_ovrlp(meta.m_ovrlp), //new edits
								 m_max_block_read(meta.m_max_block_read),
								 m_max_block_count(meta.m_max_block_count),		//edits by Samyak
								 m_tstamp_per_block(meta.m_tstamp_per_block),
                               					 m_block_ovrlp_band3(meta.m_block_ovrlp_band3),
                               					 m_block_ovrlp_band4(meta.m_block_ovrlp_band4),
                               					 m_block_ovrlp_band5(meta.m_block_ovrlp_band5),
                               					 m_enable_flag(meta.m_enable_flag), 
                               					 m_free_memory(meta.m_free_memory) { 			//edits by Samyak end
    
    }

    /** \brief Constructor for aa_filterbank_metadata that sets all values on construction. */
    aa_filterbank_metadata(const int &telescope_id,
                           const int &machine_id,
                           const int &data_type,
                           const std::string &rawdatafile,
                           const std::string &source_name,
                           const int &barycentric,
                           const int &pulsarcentric,
                           const double &az_start,
                           const double &za_start,
                           const double &src_raj,
                           const double &src_dej,
                           const double &tstart,
                           const double &tsamp,
                           const int &nbits,
                           const int &nsamples,
                           const double &fch1,
                           const double &foff,
			   const double &sign_bw, //new edits
                           const char &FREQUENCY_START,
                           const double &fchannel,
                           const char &FREQUENCY_END,
                           const int &nchans,
                           const int &nifs,
                           const double &refdm,
                           const double &period,
                           const int &strat, //new edits
                           const float &ovrlp/*, //new edits
			   const int &max_block_read,		//edits by Samyak
			   const int &max_block_count,
			   const int &tstamp_per_block,
			   const int &block_ovrlp_band3,
			   const int &block_ovrlp_band4,
			   const int &block_ovrlp_band5,
			   const int &enable_flag,
			   const int &free_memory*/)		//edits by Samyak end
      : m_az_start(az_start),
	m_za_start(za_start),
	m_src_raj(src_raj),
	m_src_dej(src_dej),
	m_tstart(tstart),
	m_tsamp(tsamp),
	m_refdm(refdm),
	m_period(period),
	m_fch1(fch1),
	m_foff(foff),
        m_sign_bw(sign_bw), //new edits
	m_fchannel(fchannel),
	m_telescope_id(telescope_id),
	m_machine_id(machine_id),
	m_data_type(data_type),
	m_barycentric(barycentric),
	m_pulsarcentric(pulsarcentric),
	m_nbits(nbits),
	m_nsamples(nsamples),
	m_nchans(nchans),
	m_nifs(nifs),
	m_FREQUENCY_START(FREQUENCY_START),
	m_FREQUENCY_END(FREQUENCY_END),
	m_rawdatafile(rawdatafile),
	m_source_name(source_name), 
	m_strat(strat), //new edits
	m_ovrlp(ovrlp)/*, 
	m_max_block_read(max_block_read),
	m_max_block_count(max_block_count),             //edits by Samyak
        m_tstamp_per_block(tstamp_per_block),
        m_block_ovrlp_band3(block_ovrlp_band3),
        m_block_ovrlp_band4(block_ovrlp_band4),
        m_block_ovrlp_band5(block_ovrlp_band5),
        m_enable_flag(enable_flag),
        m_free_memory(free_memory)*/ {                    //edits by Samyak end


 //new edits
        
    }

    /** \brief Constructor for aa_filterbank_metadata with a reduced parameter set. */
    aa_filterbank_metadata(const double &tstart,
			   const double &tsamp,
			   const int &nbits,
			   const int &nsamples,
			   const double &fch1,
			   const double &foff,
			   const double &sign_bw, //new edits
			   const int &nchans/*,
			   const int &max_block_read,           //edits by Samyak
                           const int &max_block_count,
                           const int &tstamp_per_block,
                           const int &block_ovrlp_band3,
                           const int &block_ovrlp_band4,
                           const int &block_ovrlp_band5,
                           const int &enable_flag,
                           const int &free_memory*/) : m_az_start(0),	//edits by Samyak end
						m_za_start(0),
						m_src_raj(0),
						m_src_dej(0),
						m_tstart(tstart),
						m_tsamp(tsamp),
						m_refdm(0),
						m_period(0),
						m_fch1(fch1),
						m_foff(foff),
   						m_sign_bw(sign_bw), //new edits
						m_fchannel(0),
						m_telescope_id(0),
						m_machine_id(0),
						m_data_type(0),
						m_barycentric(0),
						m_pulsarcentric(0),
						m_nbits(nbits),
						m_nsamples(nsamples),
						m_nchans(nchans),
						m_nifs(0),
						m_FREQUENCY_START(0),
						m_FREQUENCY_END(0),
						m_rawdatafile(""),
						m_source_name(""),
						m_strat(0), //new edits
						m_ovrlp(0),
						m_max_block_read(0),		//edits by Samyak
						m_max_block_count(0),
                               			m_tstamp_per_block(0),
                               			m_block_ovrlp_band3(0),
                               			m_block_ovrlp_band4(0),
                               			m_block_ovrlp_band5(0),
                               			m_enable_flag(0),
                               			m_free_memory(0) { //new edits	//edits by samyak end
     						 
    }

//brief constructor for select_data with reduced parameters
aa_filterbank_metadata(/*const double &tstart,
                           const double &tsamp,
                           const int &nbits,
                           const int &nsamples,
                           const double &fch1,
                           const double &foff,
                           const double &sign_bw, //new edits
                           const int &nchans,*/
                           const int &max_block_read,           //edits by Samyak
                           const int &max_block_count,
                           const int &tstamp_per_block,
                           const int &block_ovrlp_band3,
                           const int &block_ovrlp_band4,
                           const int &block_ovrlp_band5,
                           const int &enable_flag, 
                           const int &free_memory): /*m_az_start(0),   //edits by Samyak end
                                                m_za_start(0),
                                                m_src_raj(0),
                                                m_src_dej(0),
                                                m_tstart(tstart),
                                                m_tsamp(tsamp),
                                                m_refdm(0),
                                                m_period(0),
                                                m_fch1(fch1),
                                                m_foff(foff),
                                                m_sign_bw(sign_bw), //new edits
                                                m_fchannel(0),
                                                m_telescope_id(0),
                                                m_machine_id(0),
                                                m_data_type(0),
                                                m_barycentric(0),
                                                m_pulsarcentric(0),
                                                m_nbits(nbits),
                                                m_nsamples(nsamples),
                                                m_nchans(nchans),
                                                m_nifs(0),
                                                m_FREQUENCY_START(0),
                                                m_FREQUENCY_END(0),
                                                m_rawdatafile(""),
                                                m_source_name(""),
                                                m_strat(0), //new edits
                                                m_ovrlp(0),*/
                                                m_max_block_read(max_block_read),            //edits by Samyak
                                                m_max_block_count(max_block_count),
                                                m_tstamp_per_block(tstamp_per_block),
                                                m_block_ovrlp_band3(block_ovrlp_band3),
                                                m_block_ovrlp_band4(block_ovrlp_band4),
                                                m_block_ovrlp_band5(block_ovrlp_band5),
                                                m_enable_flag(enable_flag),
                                                m_free_memory(free_memory) { //new edits        //edits by samyak end

    }

    
    /** Destructor for aa_filterbank_metadata */
    ~aa_filterbank_metadata() {
        
    }
    
    /** \returns az_start. */
    double az_start() const {
      return m_az_start;
    }
    
    /** \returns za_start. */
    double za_start() const {
      return m_za_start;
    }
    
    /** \returns src_raj. */
    double src_raj() const {
      return m_src_raj;
    }
    
    /** \returns src_dej. */
    double src_dej() const {
      return m_src_dej;
    }
    
    /** \returns tstart. */
    double tstart() const {
      return m_tstart;
    }
    
    /** \returns tsamp. */
    double tsamp() const {
      return m_tsamp;
    }
    
    /** \returns refdm. */
    double refdm() const {
      return m_refdm;
    }
    
    /** \returns period. */
    double period() const {
      return m_period;
    }
    
    /** \returns fch1. */
    double fch1() const {
      return m_fch1;
    }
    
    /** \returns foff. */
    double foff() const {
      return m_foff;
    }    
//new edits(1)    
    /** \returns sign_bw. */
    double sign_bw() const {
      return m_sign_bw;
    }
//end(1)
    /** \returns fchannel. */
    double fchannel() const {
      return m_fchannel;
    }
    
    /** \returns telescope_id. */
    int telescope_id() const {
      return m_telescope_id;
    }
    
    /** \returns machine_id. */
    int machine_id() const {
      return m_machine_id;
    }
    
    /** \returns data_type. */
    int data_type() const {
      return m_data_type;
    }
    
    /** \returns barycentric. */
    int barycentric() const {
      return m_barycentric;
    }
    
    /** \returns pulsarcentric. */
    int pulsarcentric() const {
      return m_pulsarcentric;
    }
    
    /** \returns nbits. */
    int nbits() const {
      return m_nbits;
    }
    
    /** \returns nsamples. */
    int nsamples() const {
      return m_nsamples;
    }
    
    /** \returns nchans. */
    int nchans() const {
      return m_nchans;
    }
    
    /** \returns nifs. */
    int nifs() const {
      return m_nifs;
    }
    
    /** \returns FREQUENCY_START. */
    char FREQUENCY_START() const {
      return m_FREQUENCY_START;
    }
    
    /** \returns FREQUENCY_END. */
    char FREQUENCY_END() const {
      return m_FREQUENCY_END;
    }
    
    /** \returns rawdatafile. */
    const std::string rawdatafile() const {
      return m_rawdatafile;
    }
    
    /** \returns source_name. */
    const std::string source_name() const {
      return m_source_name;
    }
    
    /** \returns N. */
    int N() const {
      return (m_nifs * m_nchans * m_nsamples);
    }
    
    /** \returns array_index for a given samples_idx, channel_idx, and frequency_channel. */
    int array_index(const int &sample_idx, const int &channel_idx, const int &frequency_channel) const {
      return ((sample_idx * m_nifs * m_nchans)
	      + (channel_idx * m_nchans + frequency_channel));
    }
    
    /** \returns sky_frequency. */
    double sky_frequency(const int &channel_idx) const {
      return (m_fch1 + channel_idx * m_foff);
    }
//new edits(2)
    /** \returns mstrat. */
    int strat() const {
      return m_strat;
    }

    /** \returns mstrat. */
    float ovrlp() const {
      return m_ovrlp;
    }
//end(2)
//edits by Samyak 
// \returns max_block_read. 
    int max_block_read() const {
      return m_max_block_read;
    }
// \returns max_block_count. 
    int max_block_count() const {
      return m_max_block_count;
    }
// \returns tstamp_per_block. 
    int tstamp_per_block() const {
      return m_tstamp_per_block;
    }
// \returns block_ovrlp_band3. 
    int block_ovrlp_band3() const {
      return m_block_ovrlp_band3;
    }
// \returns block_ovrlp_band4. 
    int block_ovrlp_band4() const {
      return m_block_ovrlp_band4;
    }
// \returns block_ovrlp_band5. 
    int block_ovrlp_band5() const {
      return m_block_ovrlp_band5;
    }
// \returns enable_flag. 
    int enable_flag() const {
      return m_enable_flag;
    }
// \returns free_memory. 
    int free_memory() const {
      return m_free_memory;
    }

//edits by Samyak end
	

    static bool print_info(const aa_filterbank_metadata &metadata) {
      LOG(log_level::dev_debug, "FILTERBANK METADATA INFORMATION:");
      LOG(log_level::dev_debug, "az_start:\t\t" + std::to_string(metadata.az_start()));
      LOG(log_level::dev_debug, "za_start:\t\t" + std::to_string(metadata.za_start()));
      LOG(log_level::dev_debug, "src_raj:\t\t\t" + std::to_string(metadata.src_raj()));
      LOG(log_level::dev_debug, "src_dej:\t\t\t" + std::to_string(metadata.src_dej()));
      LOG(log_level::dev_debug, "tstart:\t\t\t" + std::to_string(metadata.tstart()));
      LOG(log_level::dev_debug, "tsamp:\t\t\t" + std::to_string(metadata.tsamp()));
      LOG(log_level::dev_debug, "refdm:\t\t\t" + std::to_string(metadata.refdm()));
      LOG(log_level::dev_debug, "period:\t\t\t" + std::to_string(metadata.period()));
      LOG(log_level::dev_debug, "fch1:\t\t\t" + std::to_string(metadata.fch1()));
      LOG(log_level::dev_debug, "foff:\t\t\t" + std::to_string(metadata.foff()));
      LOG(log_level::dev_debug, "sign_bw:\t\t\t" + std::to_string(metadata.sign_bw())); //new edits
      LOG(log_level::dev_debug, "fchannel:\t\t" + std::to_string(metadata.fchannel()));
      LOG(log_level::dev_debug, "telescope_id:\t\t" + std::to_string(metadata.telescope_id()));
      LOG(log_level::dev_debug, "machine_id:\t\t" + std::to_string(metadata.machine_id()));
      LOG(log_level::dev_debug, "data_type:\t\t" + std::to_string(metadata.data_type()));
      LOG(log_level::dev_debug, "barycentric:\t\t" + std::to_string(metadata.barycentric()));
      LOG(log_level::dev_debug, "pulsarcentric:\t\t" + std::to_string(metadata.pulsarcentric()));
      LOG(log_level::dev_debug, "nbits:\t\t\t" + std::to_string(metadata.nbits()));
      LOG(log_level::dev_debug, "nsamples:\t\t" + std::to_string(metadata.nsamples()));
      LOG(log_level::dev_debug, "nchans:\t\t\t" + std::to_string(metadata.nchans()));
      LOG(log_level::dev_debug, "nifs:\t\t\t" + std::to_string(metadata.nifs()));
      LOG(log_level::dev_debug, "FREQUENCY_START:\t\t" + std::to_string(metadata.FREQUENCY_START()));
      LOG(log_level::dev_debug, "FREQUENCY_END:\t\t" + std::to_string(metadata.FREQUENCY_END()));
      LOG(log_level::dev_debug, "rawdatafile:\t\t" + metadata.rawdatafile());
      LOG(log_level::dev_debug, "source_name:\t\t" + metadata.source_name());
      LOG(log_level::dev_debug, "N:\t\t\t" + std::to_string(metadata.N()));
      LOG(log_level::dev_debug, "m_strat:\t\t\t" + std::to_string(metadata.strat())); //new edits
      LOG(log_level::dev_debug, "m_ovrlp:\t\t\t" + std::to_string(metadata.ovrlp())); //new edits
      LOG(log_level::dev_debug, "max_block_read:\t\t\t" + std::to_string(metadata.max_block_read()));		//edits by Samyak
      LOG(log_level::dev_debug, "max_block_count:\t\t\t" + std::to_string(metadata.max_block_count()));
      LOG(log_level::dev_debug, "tstamp_per_block:\t\t\t" + std::to_string(metadata.tstamp_per_block()));
      LOG(log_level::dev_debug, "block_ovrlp_band3:\t\t\t" + std::to_string(metadata.block_ovrlp_band3()));
      LOG(log_level::dev_debug, "block_ovrlp_band4:\t\t\t" + std::to_string(metadata.block_ovrlp_band4()));
      LOG(log_level::dev_debug, "block_ovrlp_band5:\t\t\t" + std::to_string(metadata.block_ovrlp_band5()));
      LOG(log_level::dev_debug, "enable_flag:\t\t\t" + std::to_string(metadata.enable_flag()));
      LOG(log_level::dev_debug, "free_memory:\t\t\t" + std::to_string(metadata.free_memory())); 		// edits by Samyak end

      return true;
    }
    
  private:
    /** For data descriptions, please see SIGPROC interface control document as referenced in the class. */
    double m_az_start;
    double m_za_start;
    double m_src_raj;
    double m_src_dej;
    double m_tstart;
    double m_tsamp;
    double m_refdm;
    double m_period;
    double m_fch1;
    double m_foff;
    double m_sign_bw; //new edits
    double m_fchannel;
    
    int m_telescope_id;
    int m_machine_id;
    int m_data_type;
    int m_barycentric;
    int m_pulsarcentric;
    int m_nbits;
    int m_nsamples;
    int m_nchans;
    int m_nifs;
    int m_strat; //new edits
    float m_ovrlp; //new edits
    int m_max_block_read;		//edits by Samyak
    int m_max_block_count;
    int m_tstamp_per_block;
    int m_block_ovrlp_band3;
     int m_block_ovrlp_band4;
    int m_block_ovrlp_band5;
    int m_enable_flag;
    int m_free_memory;
    char m_FREQUENCY_START;
    char m_FREQUENCY_END;	//edits by Samyak end
    
    std::string m_rawdatafile;
    std::string m_source_name;
  };
} // namespace astroaccelerate
  
#endif // ASTRO_ACCELERATE_AA_FILTERBANK_METADATA_HPP
