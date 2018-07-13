/**
 * \file
 * \date  14 June 2018
 * \brief Open a possible compressed file for reading
 */
#ifndef _OpenFile_h
#define _OpenFile_h
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/copy.hpp>

#ifdef USE_ZLIB
#include <boost/iostreams/filter/gzip.hpp>
#endif

#ifdef USE_BZIP2
#include <boost/iostreams/filter/bzip2.hpp>
#endif

#ifdef USE_LZMA
#include <boost/iostreams/filter/lzma.hpp>
#endif

/**
 * Opens a file for reading
 *
 * Determines the type of file and sets up the appropriate
 * boost::iostreams::filtering_stream.
 *
 * \param[in]      filename Name of file
 * \param[in,out]  in       boost filtering_stream
 * \param[in,out]  ifp      ifstream
 *
 * \return \c TRUE if the file could be open, \c FALSE otherwise
 */
bool openfile(const char *filename,
              boost::iostreams::filtering_stream<boost::iostreams::input> *in,
              std::ifstream *ifp);
#endif
