/*
 * (C) Copyright 1996-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

/// @author Willem Deconinck
/// @date   Nov 2013

#ifndef HaloExchange_h
#define HaloExchange_h

#include <string>
#include <vector>
#include <stdexcept>
#include "eckit/memory/SharedPtr.h"
#include "eckit/memory/Owned.h"
#include "eckit/mpi/Exceptions.h"
#include "eckit/exception/Exceptions.h"
#include "atlas/internals/Debug.h"
#include "atlas/array/ArrayView.h"
#include "atlas/parallel/mpi/mpi.h"

namespace atlas {
namespace parallel {

class HaloExchange: public eckit::Owned {

public: // types
    typedef eckit::SharedPtr<HaloExchange> Ptr;
public:
  HaloExchange();
  HaloExchange(const std::string& name);
  virtual ~HaloExchange() {}

public: // methods

  const std::string& name() const { return name_; }

  void setup( const int part[],
              const int remote_idx[], const int base,
              size_t size );

  template <typename DATA_TYPE>
  void execute( DATA_TYPE field[], const size_t var_strides[], const size_t var_shape[], size_t var_rank ) const;

  template <typename DATA_TYPE>
  void execute( DATA_TYPE field[], size_t nb_vars ) const;

  template <typename DATA_TYPE, int RANK>
  void execute( array::ArrayView<DATA_TYPE,RANK>& field ) const;

private: // methods

  void create_mappings( std::vector<int>& send_map, std::vector<int>& recv_map, size_t nb_vars ) const;

  template<int N, int P>
  void create_mappings_impl( std::vector<int>& send_map, std::vector<int>& recv_map, size_t nb_vars ) const;

  size_t index(size_t i, size_t j, size_t k, size_t ni, size_t nj, size_t nk) const { return( i + ni*( j + nj*k) ); }

  size_t index(size_t i, size_t j, size_t ni, size_t nj) const { return( i + ni*j ); }


  template< typename DATA_TYPE>
  void pack_send_buffer( const DATA_TYPE field[],
                         const size_t var_strides[],
                         const size_t var_shape[],
                         size_t var_rank,
                         DATA_TYPE send_buffer[] ) const;

  template< typename DATA_TYPE>
  void unpack_recv_buffer(const DATA_TYPE recv_buffer[],
                          DATA_TYPE field[],
                          const size_t var_strides[],
                          const size_t var_shape[],
                          size_t var_rank ) const;

  template<typename DATA_TYPE, int RANK>
  void var_info( const array::ArrayView<DATA_TYPE,RANK>& arr,
                 std::vector<size_t>& varstrides,
                 std::vector<size_t>& varshape ) const;

private: // data
  std::string name_;
  bool is_setup_;

  int               sendcnt_;
  int               recvcnt_;
  std::vector<int>  sendcounts_;
  std::vector<int>  senddispls_;
  std::vector<int>  recvcounts_;
  std::vector<int>  recvdispls_;
  std::vector<int>  sendmap_;
  std::vector<int>  recvmap_;
  int parsize_;

  int nproc;
  int myproc;

  std::vector<int> bounds_;
  int par_bound_;

};


template<typename DATA_TYPE>
void HaloExchange::execute(DATA_TYPE field[], const size_t var_strides[], const size_t var_shape[], size_t var_rank ) const
{
  if( ! is_setup_ )
  {
    throw eckit::SeriousBug("HaloExchange was not setup",Here());
  }

  int tag=1;
  size_t var_size = std::accumulate(var_shape,var_shape+var_rank,1,std::multiplies<size_t>());
  int send_size  = sendcnt_ * var_size;
  int recv_size  = recvcnt_ * var_size;

  std::vector<DATA_TYPE  > send_buffer(send_size);
  std::vector<DATA_TYPE  > recv_buffer(recv_size);
  std::vector<int        > send_displs(nproc    );
  std::vector<int        > recv_displs(nproc    );
  std::vector<int        > send_counts(nproc    );
  std::vector<int        > recv_counts(nproc    );
  std::vector<MPI_Request> send_req   (nproc    );
  std::vector<MPI_Request> recv_req   (nproc    );

  for (size_t jproc=0; jproc<nproc; ++jproc)
  {
    send_counts[jproc] = sendcounts_[jproc]*var_size;
    recv_counts[jproc] = recvcounts_[jproc]*var_size;
    send_displs[jproc] = senddispls_[jproc]*var_size;
    recv_displs[jproc] = recvdispls_[jproc]*var_size;
  }


  /// Let MPI know what we like to receive
  for( size_t jproc=0; jproc<nproc; ++jproc )
  {
    if(recv_counts[jproc] > 0)
    {
      ECKIT_MPI_CHECK_RESULT( MPI_Irecv( &recv_buffer[recv_displs[jproc]] , recv_counts[jproc],
        eckit::mpi::datatype<DATA_TYPE>(), jproc, tag, eckit::mpi::comm(), &recv_req[jproc] ) );
    }
  }

  /// Pack
  pack_send_buffer(field,var_strides,var_shape,var_rank,send_buffer.data());

  /// Send
  for( size_t jproc=0; jproc<nproc; ++jproc )
  {
    if(send_counts[jproc] > 0)
    {
      ECKIT_MPI_CHECK_RESULT( MPI_Isend( &send_buffer[send_displs[jproc]], send_counts[jproc],
        eckit::mpi::datatype<DATA_TYPE>(), jproc, tag, eckit::mpi::comm(), &send_req[jproc] ) );
    }
  }

  /// Wait for receiving to finish
  for (size_t jproc=0; jproc<nproc; ++jproc)
  {
    if( recvcounts_[jproc] > 0)
    {
      ECKIT_MPI_CHECK_RESULT( MPI_Wait(&recv_req[jproc], MPI_STATUS_IGNORE ) );
    }
  }

  /// Unpack
  unpack_recv_buffer(recv_buffer.data(),field,var_strides,var_shape,var_rank);

  /// Wait for sending to finish
  for (size_t jproc=0; jproc<nproc; ++jproc)
  {
    if( sendcounts_[jproc] > 0)
    {
      ECKIT_MPI_CHECK_RESULT( MPI_Wait(&send_req[jproc], MPI_STATUS_IGNORE ) );
    }
  }

}




template<typename DATA_TYPE>
void HaloExchange::pack_send_buffer( const DATA_TYPE field[],
                                     const size_t var_strides[],
                                     const size_t var_shape[],
                                     size_t var_rank,
                                     DATA_TYPE send_buffer[] ) const
{
  size_t ibuf = 0;
  size_t send_stride = var_strides[0]*var_shape[0];

  switch( var_rank )
  {
  case 1:
    for( size_t p=0; p<sendcnt_; ++p)
    {
      const size_t pp = send_stride*sendmap_[p];
      for( size_t i=0; i<var_shape[0]; ++i )
        send_buffer[ibuf++] = field[pp+i*var_strides[0]];
    }
    break;
  case 2:
    for( size_t p=0; p<sendcnt_; ++p)
    {
      const size_t pp = send_stride*sendmap_[p];
      for( size_t i=0; i<var_shape[0]; ++i )
      {
        for( size_t j=0; j<var_shape[1]; ++j )
        {
          send_buffer[ibuf++] = field[pp+i*var_strides[0]+j*var_strides[1]];
        }
      }
    }
    break;
  case 3:
    for( size_t p=0; p<sendcnt_; ++p)
    {
      const size_t pp = send_stride*sendmap_[p];
      for( size_t i=0; i<var_shape[0]; ++i )
      {
        for( size_t j=0; j<var_shape[1]; ++j )
        {
          for( size_t k=0; k<var_shape[2]; ++k )
          {
            send_buffer[ibuf++] =
              field[ pp+i*var_strides[0]+j*var_strides[1]+k*var_strides[2]];
          }
        }
      }
    }
    break;
  case 4:
    for( size_t p=0; p<sendcnt_; ++p)
    {
      const size_t pp = send_stride*sendmap_[p];
      for( size_t i=0; i<var_shape[0]; ++i )
      {
        for( size_t j=0; j<var_shape[1]; ++j )
        {
          for( size_t k=0; k<var_shape[2]; ++k )
          {
           for( size_t l=0; l<var_shape[3]; ++l )
           {
            send_buffer[ibuf++] =
              field[ pp+i*var_strides[0]+j*var_strides[1]+k*var_strides[2]+l*var_strides[3]];
           }
          }
        }
      }
    }
    break;
  default:
    NOTIMP;
  }
}

template<typename DATA_TYPE>
void HaloExchange::unpack_recv_buffer( const DATA_TYPE recv_buffer[],
                                       DATA_TYPE field[],
                                       const size_t var_strides[],
                                       const size_t var_shape[],
                                       size_t var_rank ) const
{
//  bool field_changed = false;
//  DATA_TYPE tmp;
  size_t ibuf = 0;
  size_t recv_stride = var_strides[0]*var_shape[0];

  switch( var_rank )
  {
  case 1:
    for( size_t p=0; p<recvcnt_; ++p)
    {
      const size_t pp = recv_stride*recvmap_[p];
      for( size_t i=0; i<var_shape[0]; ++i)
      {
//        tmp = field[ pp + i*var_strides[0] ];
        field[ pp + i*var_strides[0] ] = recv_buffer[ibuf++];
//        if( tmp != field[ pp + i*var_strides[0] ] )
//          field_changed = true;
      }
    }
    break;
  case 2:
    for( size_t p=0; p<recvcnt_; ++p)
    {
      const size_t pp = recv_stride*recvmap_[p];
      for( size_t i=0; i<var_shape[0]; ++i )
      {
        for( size_t j=0; j<var_shape[1]; ++j )
        {
//          tmp = field[ pp + i*var_strides[0] + j*var_strides[1] ];
          field[ pp + i*var_strides[0] + j*var_strides[1] ]
              = recv_buffer[ibuf++];
//          if( field[ pp + i*var_strides[0] + j*var_strides[1] ] != tmp )
//            field_changed = true;
        }
      }
    }
    break;
  case 3:
    for( size_t p=0; p<recvcnt_; ++p)
    {
      const size_t pp = recv_stride*recvmap_[p];
      for( size_t i=0; i<var_shape[0]; ++i )
      {
        for( size_t j=0; j<var_shape[1]; ++j )
        {
          for( size_t k=0; k<var_shape[2]; ++k )
          {
//            tmp = field[ pp + i*var_strides[0] + j*var_strides[1] + k*var_strides[2] ];
            field[ pp + i*var_strides[0] + j*var_strides[1] + k*var_strides[2] ]
                = recv_buffer[ibuf++];
//            if( field[ pp + i*var_strides[0] + j*var_strides[1] + k*var_strides[2] ] != tmp )
//              field_changed = true;
          }
        }
      }
    }
    break;
  case 4:
    for( size_t p=0; p<recvcnt_; ++p)
    {
      const size_t pp = recv_stride*recvmap_[p];
      for( size_t i=0; i<var_shape[0]; ++i )
      {
        for( size_t j=0; j<var_shape[1]; ++j )
        {
          for( size_t k=0; k<var_shape[2]; ++k )
          {
           for( size_t l=0; l<var_shape[3]; ++l )
           {
//            tmp = field[ pp + i*var_strides[0] + j*var_strides[1] + k*var_strides[2] +l*var_strides[3] ];
            field[ pp + i*var_strides[0] + j*var_strides[1] + k*var_strides[2] + l*var_strides[3] ]
                = recv_buffer[ibuf++];
//            if( field[ pp + i*var_strides[0] + j*var_strides[1] + k*var_strides[2] +l*var_strides[3] ] != tmp )
//              field_changed = true;
           }
          }
        }
      }
    }
    break;
  default:
    NOTIMP;
  }
//  if( !field_changed )
//    std::cout << "WARNING: halo-exchange did not change field" << std::endl;
}

template<typename DATA_TYPE>
void HaloExchange::execute( DATA_TYPE field[], size_t nb_vars ) const
{
  size_t strides[] = {1};
  size_t shape[] = {nb_vars};
  execute( field, strides, shape, 1);
}


template<typename DATA_TYPE, int RANK>
void HaloExchange::var_info( const array::ArrayView<DATA_TYPE,RANK>& arr,
                             std::vector<size_t>& varstrides,
                             std::vector<size_t>& varshape ) const
{
  size_t rank = std::max(1,RANK-1) ;
  varstrides.resize(rank);
  varshape.resize(rank);
  if( RANK>1 )
  {
    varstrides.assign(arr.strides()+1,arr.strides()+RANK);
    varshape.assign(arr.shape()+1,arr.shape()+RANK);
  }
  else
  {
    varstrides[0] = arr.strides()[0];
    varshape[0] = 1;
  }
}

template <typename DATA_TYPE, int RANK>
void HaloExchange::execute( array::ArrayView<DATA_TYPE,RANK>& field ) const
{
  //if( field.shape(0) == parsize_)
  if( true ){
    std::vector<size_t> varstrides, varshape;
    var_info( field, varstrides, varshape );
    execute( field.data(), varstrides.data(), varshape.data(), varstrides.size() );
  }
  else
  {
    Log::error() << "Passed field with rank " << RANK << " and shape(0) " << field.shape(0) << std::endl;
    Log::error() << "parsize_ = " << parsize_ << std::endl;
    NOTIMP; // Need to implement with parallel ranks > 1
  }
}


// ------------------------------------------------------------------
// C wrapper interfaces to C++ routines
extern "C"
{
  HaloExchange* atlas__HaloExchange__new ();
  void atlas__HaloExchange__delete (HaloExchange* This);
  void atlas__HaloExchange__setup (HaloExchange* This, int part[], int remote_idx[], int base, int size);
  void atlas__HaloExchange__execute_strided_int (HaloExchange* This, int field[], int var_strides[], int var_shape[], int var_rank);
  void atlas__HaloExchange__execute_strided_long (HaloExchange* This, long field[], int var_strides[], int var_shape[], int var_rank);
  void atlas__HaloExchange__execute_strided_float (HaloExchange* This, float field[], int var_strides[], int var_shape[], int var_rank);
  void atlas__HaloExchange__execute_strided_double (HaloExchange* This, double field[], int var_strides[], int var_shape[], int var_rank);
  void atlas__HaloExchange__execute_int (HaloExchange* This, int field[], int var_rank);
  void atlas__HaloExchange__execute_float (HaloExchange* This, float field[], int var_rank);
  void atlas__HaloExchange__execute_double (HaloExchange* This, double field[], int var_rank);

}
// ------------------------------------------------------------------

} // namespace parallel
} // namespace atlas

#endif // HaloExchange_h