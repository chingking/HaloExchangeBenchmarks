
/*
Copyright (c) 2012, MAURO BIANCO, UGO VARETTO, SWISS NATIONAL SUPERCOMPUTING CENTRE (CSCS)
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.
    * Neither the name of the Swiss National Supercomputing Centre (CSCS) nor the
      names of its contributors may be used to endorse or promote products
      derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL MAURO BIANCO, UGO VARETTO, OR 
SWISS NATIONAL SUPERCOMPUTING CENTRE (CSCS), BE LIABLE FOR ANY
DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/


#ifndef _HALO_EXCHANGE_H_
#define _HALO_EXCHANGE_H_

#include <Halo_Exchange_2D.h>
#include <Halo_Exchange_2D_DT.h>
#include <Halo_Exchange_3D.h>
#include <Halo_Exchange_3D_DT.h>
#include <proc_grids_2D.h>
#include <proc_grids_3D.h>
#include <utils/boollist.h>

#ifdef MANUAL_PACKING
#error("Manual Packing is now turned on by setting versions to GCL::version_manual (or, equivalently) 2")
#endif

//#else
#include <descriptors_dt.h>
#include <descriptors_dt_whole.h>
#include <descriptors.h>
#include <descriptors_manual_gpu.h>
#include <descriptor_generic_manual.h>
 //#endif

#include <translate.h>
#include <field_on_the_fly.h>

namespace GCL {

  namespace _impl {
    template <int D, typename GT, int version>
    struct get_pattern;

    template <typename GT>
    struct get_pattern<2,GT,0> {
      typedef Halo_Exchange_2D<GT> type;
    };

    template <typename GT>
    struct get_pattern<2,GT,1> {
      typedef Halo_Exchange_2D_DT<GT> type;
    };

    template <typename GT>
    struct get_pattern<2,GT,2> {
      typedef Halo_Exchange_2D<GT> type;
    };

    template <typename GT>
    struct get_pattern<3,GT,0> {
      typedef Halo_Exchange_3D<GT> type;
    };

    template <typename GT>
    struct get_pattern<3,GT,1> {
      typedef Halo_Exchange_3D_DT<GT> type;
    };

    template <typename GT>
    struct get_pattern<3,GT,2> {
      typedef Halo_Exchange_3D<GT> type;
    };

    template <int D>
    struct get_grid;

    template <>
    struct get_grid<2> {
      typedef MPI_2D_process_grid_t<gcl_utils::boollist<2> > type; 
    };

    template <>
    struct get_grid<3> {
      typedef MPI_3D_process_grid_t<gcl_utils::boollist<3> > type;
    };
    
  }


  /** 
     This is the main class for the halo exchange pattern in the case
     in which the data pointers are not known before hand and they are
     passed to the pattern when packing and unpacking is needed.

     The interface requires two layout maps ( \link GCL::layout_map
     \endlink ) one for to specify the data layout, the other to
     specify the relation between data layout and processor grid. This
     is an important asepct that will be explained here and also in
     the introduction.

     The First layout map to be passed to the pattern class is the
     data layout. The user defines a convention in which the
     dimensions of the data fields are ordered logically depending on
     the application and/or user preferences. For instance, we can
     call the dimensions in this application order i, j, and k. The
     the layout map in this case specifies in what position each
     dimension is in the increasing stride order. For instance:

     \code
     GCL::layout_map<1,0,2>
     \endcode

     Indicates that the first dimension in the data (i) is the second
     in the increasing stride order, while the second (j) is actually
     the first one (has stride 1). Finally the k dimension is the one
     with the highest stride and is in the third position (2). The
     layout in memory of the user data can be specified then as 'jik'.

     Similarly, the second template argument in the halo exchange
     pattern is the map between data coordinates and the processor
     grid coordinates. The following layout specification

     \code
     GCL::layout_map<1,0,2>
     \endcode

     would mean: The first dimension in data matches with the second
     dimension of the computing grid, the second is the first, and the
     third one is the third one.

     Let's consider a 2D case at first, to show an additional
     example. Suppose user data is thought to be ij, meaning the user
     think to i as the first coordinate and j to the
     second. Alternatively the user would use (i,j) to indicate an
     element in the data array. The layout is C like, so that j is
     actuallly the first coordinate in the increasing stride ordering,
     and i the second.

     The first template argument to the pattern would then be
     \code
     GCL::layout_map<1,0>
     \endcode

     The second template argument is still a \link GCL::layout_map
     \endlink , but this time it indicates the mapping between data
     and processor grid. The data is still condidered in the user
     convention.

     Suppose the processor gris has size PIxPJ, and these may be the
     MPI communicator sizes (coords[0]=PI, coords[1]=PJ). Now, we want
     to say that the first dimension on data (first in the user
     convention, not int the increasing stride order) 'extends' to the
     computing gris, or that the first dimension in the data
     correspons to the first dimension in the computing grid. Let's
     consider a 2x1 proc grid, and the first dimension of the data
     being the rows (i) and the second the column (j). In this case we
     are thinking to a distribution like this:

     \code
        >j>>
       ------
     v |0123|
     i |1234|  Proc 0,0
     v |2345|
     v |3456|
       ------

        >j>>
       ------
     v |4567|
     i |5678|  Proc 1,0
     v |6789|
     v |7890|
       ------
     \endcode

     In this case the map between data and the processor grid is:
     \code
     GCL::layout_map<0,1>
     \endcode

     On the other hand, having specified 
     \code
     GCL::layout_map<1,0>
     \endcode
     for this map, would imply a layout/distribution like the following:

     \code
        >j>>                 >j>>
       ------               ------
     v |0123|             v |4567|
     i |1234|  Proc 0,0;  i |5678|  Proc 1,0
     v |2345|             v |6789|
     v |3456|             v |7890|
       ------               ------
     \endcode

     Where the second dimension in the data correspond to the fist
     dimension in the processor grid. Again, the data coordinates
     ordering is the one the user choose to be the logical order in
     the application, not the increasing stride order.

     To find an example in 2D refer to 
     test_halo_exchange_2D.cpp

     while a 3D example can be found in 
     test_halo_exchange.cpp

     The other template arguments are the type of the elements
     contained in the data arrays and the number of dimensions of
     data.

     \tparam layout_map Layout_map \link GCL::layout_map \endlink specifying the data layout as the position of each dimension of the user data in the increasing stride order
     \tparam layout2proc_map_abs Layout_map \link GCL::layout_map \endlink specifying which dimension in the data corresponds to the which dimension in the processor grid
     \tparam DataType Value type the elements int the arrays
     \tparam DIMS Number of dimensions of data arrays (equal to the dimension of the processor grid)
     \tparam GCL_ARCH Specification of the "architecture", that is the place where the data to be exchanged is. Possible coiches are defined in L3/include/gcl_arch.h .
   */
  template <typename layout_map,
            typename layout2proc_map_abs,
            typename DataType,
            int DIMS,
            typename Gcl_Arch = gcl_cpu,
            int version=0>
  class halo_exchange_dynamic_ut {

  private:
    typedef typename layout_transform<layout_map, layout2proc_map_abs>::type layout2proc_map;
    
  public:
    /**
       Type of the computin grid associated to the pattern
    */
    typedef typename _impl::get_grid<DIMS>::type grid_type;

    /**
       Type of the Level 3 pattern used. This is available only if the pattern uses a Level 3 pattern. 
       In the case the implementation is not using L3, the type is not available.
    */
    typedef typename _impl::get_pattern<DIMS, grid_type, version>::type pattern_type;

  private:
    hndlr_dynamic_ut<DataType,
                     DIMS, 
                     pattern_type,
                     layout2proc_map,
                     Gcl_Arch,
                     version> hd;

    halo_exchange_dynamic_ut(halo_exchange_dynamic_ut const &) {}

    typename grid_type::period_type periodicity;

#ifdef GCL_TRACE
    int pattern_tag;
#endif

  public:

    /** constructor that takes the periodicity (mathich the \link
        boollist_concept \endlink concept, and the MPI CART
        communicator in DIMS (specified as template argument to the
        pattern) dimensions of the processing grid. the periodicity is
        specified in the order chosen by the programmer for the data,
        as in the rest of the application. It is up tp the
        construnctor implementation to translate it into the right
        order depending on the GCL::layout_map passed to the class.

        \param[in] c Periodicity specification as in \link boollist_concept \endlink
        \param[in] comm MPI CART communicator with dimension DIMS (specified as template argument to the pattern).
    */

    template <typename layout, int DIM>
    struct proc_map {};

    template <typename layout>
    struct proc_map<layout, 2>{
        static std::vector<int> map(){
            std::vector<int> m(2);
            m[0] = layout::template at<0>();
            m[1] = layout::template at<1>();
            return m;
        }
    };

    template <typename layout>
    struct proc_map<layout, 3>{
        static std::vector<int> map(){
            std::vector<int> m(3);
            m[0] = layout::template at<0>();
            m[1] = layout::template at<1>();
            m[2] = layout::template at<2>();
            return m;
        }
    };

    explicit halo_exchange_dynamic_ut(typename grid_type::period_type const &c, MPI_Comm comm)
      : hd(c.template permute<layout2proc_map_abs>(), comm), periodicity(c)
    {
    }

    /** Function to rerturn the L3 level pattern used inside the pattern itself.

        \return The pattern al level 3 used to exchange data
    */
    pattern_type const& pattern() const {return hd.pattern();}

    /**
       Function to setup internal data structures for data exchange and preparing eventual underlying layers

       \param max_fields_n Maximum number of data fields that will be passed to the communication functions
    */
    void setup(int max_fields_n) {
      hd.setup(max_fields_n);
#ifdef GCL_TRACE
      stats_collector<DIMS>::instance()->init( hd.pattern().proc_grid().communicator );
      std::vector<int> map = proc_map<layout_map, DIMS>::map();
      int coords[DIMS];
      int dims[DIMS];
      hd.pattern().proc_grid().coords(coords[0], coords[1], coords[2]);
      hd.pattern().proc_grid().dims(dims[0], dims[1], dims[2]);
      pattern_tag = stats_collector<DIMS>::instance()->add_pattern(
                                                          Pattern<DIMS>( pt_dynamic,
                                                                         hd.halo.halos,
                                                                         map,
                                                                         periodicity,
                                                                         coords,
                                                                         dims)
                                                       );
      hd.set_pattern_tag(pattern_tag);
#endif
    }

    /**
       Function to register halos with the pattern. The registration
       happens specifing the ordiring of the dimensions as the user
       defines it. For instance in the example given at the
       introduction of the class, the ordering will be 0 for i and 1
       for j, no matter of the layout_map parameter passed during
       class instantiation.
 
       \tparam DI index of the dimension to be set relative to the logical ordering chosen in the application, not the increasing stride ordering.
       \param[in] minus Please see field_descriptor, halo_descriptor or \link MULTI_DIM_ACCESS \endlink for details
       \param[in] plus Please see field_descriptor, halo_descriptor or \link MULTI_DIM_ACCESS \endlink for details
       \param[in] begin Please see field_descriptor, halo_descriptor or \link MULTI_DIM_ACCESS \endlink for details
       \param[in] end Please see field_descriptor, halo_descriptor or \link MULTI_DIM_ACCESS \endlink for details
       \param[in] t_len Please see field_descriptor, halo_descriptor or \link MULTI_DIM_ACCESS \endlink for details
     */
    template <int DI>
    void add_halo(int minus, int plus, int begin, int end, int t_len) {
      
      hd.halo.add_halo(layout_map::template at<DI>(), minus, plus, begin, end, t_len);
    }

    /**
       Function to pack data to be sent
       
       \param[in] _fields data fields to be packed
    */
#ifdef __GXX_EXPERIMENTAL_CXX0X__
    template <typename... FIELDS>
    void pack(const FIELDS&... _fields) const {
      hd.pack(_fields... );
    }
#else
#define MACRO_IMPL(z, n, _)                                             \
    template <BOOST_PP_ENUM_PARAMS_Z(z, BOOST_PP_INC(n), typename FIELD)> \
    void pack(BOOST_PP_ENUM_BINARY_PARAMS_Z(z, BOOST_PP_INC(n), FIELD, const &_field)) const { \
      hd.pack(BOOST_PP_ENUM_PARAMS_Z(z, BOOST_PP_INC(n), _field));      \
    }
  
    BOOST_PP_REPEAT(GCL_MAX_FIELDS, MACRO_IMPL, all)
#undef MACRO_IMPL
#endif

    /**
       Function to unpack received data

       \param[in] _fields data fields where to unpack data
    */
#ifdef __GXX_EXPERIMENTAL_CXX0X__
    template <typename... FIELDS>
    void unpack(const FIELDS&... _fields) const {
      hd.unpack(_fields... );
    }
#else
#define MACRO_IMPL(z, n, _)                                             \
    template <BOOST_PP_ENUM_PARAMS_Z(z, BOOST_PP_INC(n), typename FIELD)> \
    void unpack(BOOST_PP_ENUM_BINARY_PARAMS_Z(z, BOOST_PP_INC(n), FIELD, const &_field)) const { \
      hd.unpack(BOOST_PP_ENUM_PARAMS_Z(z, BOOST_PP_INC(n), _field));    \
    }
    
    BOOST_PP_REPEAT(GCL_MAX_FIELDS, MACRO_IMPL, all)
#undef MACRO_IMPL
#endif

    /**
       Function to unpack received data

       \param[in] fields vector with data fields pointers to be packed from
    */
    void pack(std::vector<DataType*> const& fields) {
#ifdef GCL_TRACE
      double start_time = MPI_Wtime();
#endif
      hd.pack(fields);
#ifdef GCL_TRACE
#ifdef __CUDACC__
      cudaDeviceSynchronize();
#endif
      double end_time = MPI_Wtime();
      stats_collector<DIMS>::instance()->add_event(ExchangeEvent(ee_pack, start_time, end_time, fields.size(), pattern_tag));
#endif
    }

    /**
       Function to unpack received data

       \param[in] fields vector with data fields pointers to be unpacked into
    */
    void unpack(std::vector<DataType*> const& fields) {
#ifdef GCL_TRACE
      double start_time = MPI_Wtime();
#endif
      hd.unpack(fields);
#ifdef GCL_TRACE
#ifdef __CUDACC__
      cudaDeviceSynchronize();
#endif
      double end_time = MPI_Wtime();
      stats_collector<DIMS>::instance()->add_event(ExchangeEvent(ee_unpack, start_time, end_time, fields.size(), pattern_tag));
#endif
    }

    /**
       function to trigger data exchange

       Note: when the start_exchange() + wait() combination is used, the exchange() method should not be used, and vice versa.
    */
    void exchange() {
#ifdef GCL_TRACE
      double start_time = MPI_Wtime();
#endif
      hd.exchange();
#ifdef GCL_TRACE
      double end_time = MPI_Wtime();
      stats_collector<DIMS>::instance()->add_event(ExchangeEvent(ee_exchange, start_time, end_time, 0, pattern_tag));
#endif
    }

      void post_receives() {
#ifdef GCL_TRACE
      double start_time = MPI_Wtime();
#endif
          hd.post_receives();
#ifdef GCL_TRACE
      double end_time = MPI_Wtime();
      stats_collector<DIMS>::instance()->add_event(ExchangeEvent(ee_post_receives, start_time, end_time, 0, pattern_tag));
#endif
      }

      void do_sends() {
#ifdef GCL_TRACE
      double start_time = MPI_Wtime();
#endif
          hd.do_sends();
#ifdef GCL_TRACE
      double end_time = MPI_Wtime();
      stats_collector<DIMS>::instance()->add_event(ExchangeEvent(ee_do_sends, start_time, end_time, 0, pattern_tag));
#endif
      }

    /**
       function to trigger data exchange initiation when using split-phase communication.

       Note: when the start_exchange() + wait() combination is used, the exchange() method should not be used, and vice versa.
    */
    void start_exchange() {
#ifdef GCL_TRACE
      double start_time = MPI_Wtime();
#endif
      hd.start_exchange();
#ifdef GCL_TRACE
      double end_time = MPI_Wtime();
      stats_collector<DIMS>::instance()->add_event(ExchangeEvent(ee_start_exchange, start_time, end_time, 0, pattern_tag));
#endif
    }

    /**
       function to trigger data exchange

       Note: when the start_exchange() + wait() combination is used, the exchange() method should not be used, and vice versa.
    */
    void wait() {
#ifdef GCL_TRACE
      double start_time = MPI_Wtime();
#endif
      hd.wait();
#ifdef GCL_TRACE
      double end_time = MPI_Wtime();
      stats_collector<DIMS>::instance()->add_event(ExchangeEvent(ee_wait, start_time, end_time, 0, pattern_tag));
#endif
    }

  };


  template <int I>
  struct pick_version;

  template <>
  struct pick_version<2> {
    static const int value = GCL::version_mpi_pack;
  };

  template <>
  struct pick_version<3> {
    static const int value = GCL::version_manual;
  };

  /**
     This is the main class for the halo exchange pattern in the case
     in which the data pointers, data types, and shapes are not known
     before hand and they are passed to the pattern when packing and
     unpacking is needed. This is the most generic pattern available
     in GCL for halo exchange.

     The interface requires one layout map ( \link GCL::layout_map
     \endlink ) to specify the relation between data layout and
     processor grid. as follows The following layout specification

     \code
     GCL::layout_map<1,0,2>
     \endcode

     would mean: The first dimension in data matches with the second
     dimension of the computing grid, the second is the first, and the
     third one is the third one.

     \tparam layout2proc_map_abs Layout_map \link GCL::layout_map \endlink specifying which dimension in the data corresponds to the which dimension in the processor grid
     \tparam DIMS Number of dimensions of data arrays (equal to the dimension of the processor grid)
     \tparam GCL_ARCH Specification of the "architecture", that is the place where the data to be exchanged is. Possible coiches are defined in L3/include/gcl_arch.h .
   */
  template <typename layout2proc_map,
            int DIMS,
            typename Gcl_Arch = gcl_cpu,
            int version = pick_version<DIMS>::value>
  class halo_exchange_generic_base {

    //    typedef typename layout_transform<layout_map, layout2proc_map_abs>::type layout2proc_map;

  public:

    /**
       Type of the computin grid associated to the pattern
    */
    typedef typename _impl::get_grid<DIMS>::type grid_type;

    /**
       Type of the Level 3 pattern used. This is available only if the pattern uses a Level 3 pattern. 
       In the case the implementation is not using L3, the type is not available.
    */
    typedef typename _impl::get_pattern<DIMS, grid_type,version>::type pattern_type;

  private:
    hndlr_generic<DIMS, 
                  pattern_type,
                  layout2proc_map,
                  Gcl_Arch, version> hd;

  public:

    /** constructor that takes the periodicity (mathich the \link
        boollist_concept \endlink concept, and the MPI CART
        communicator in DIMS (specified as template argument to the
        pattern) dimensions of the processing grid. the periodicity is
        specified in the order chosen by the programmer for the data,
        as in the rest of the application. It is up tp the
        construnctor implementation to translate it into the right
        order depending on the GCL::layout_map passed to the class.

        \param[in] c Periodicity specification as in \link boollist_concept \endlink
        \param[in] comm MPI CART communicator with dimension DIMS (specified as template argument to the pattern).
    */
    explicit halo_exchange_generic_base(typename grid_type::period_type const &c, MPI_Comm comm)
      : hd(grid_type(c.template permute<layout2proc_map>(), comm))
    { }

    explicit halo_exchange_generic_base(grid_type const& g)
      : hd(g)
    { }

    // halo_exchange_generic(halo_exchange_generic const &src) 
    //   : 
    // {}

    /** Function to rerturn the L3 level pattern used inside the pattern itself.

        \return The pattern al level 3 used to exchange data
    */
    pattern_type const& pattern() const {return hd.pattern();}

    /**
       Function to setup internal data structures for data exchange and preparing eventual underlying layers. The sizes are computed mulplitplying the entities.

       \param max_fields_n Maximum, or sufficiently large (if compensated by other parameters), number of data fields that will be passed to the communication functions
       \param halo_example the maximum, or a suffciently large, halo from which the pattern can determine the needs of the data exchange
       \param typesize Maximum, or sufficiently large, size fo the types of values to be exchanged.
    */
    template <typename DataType, typename layomap, template <typename> class _traits>
    void setup(int max_fields_n, 
               field_on_the_fly<DataType, layomap, _traits> const & halo_example,
               int typesize) {
      hd.setup(max_fields_n, halo_example, typesize);
    }

    /**
       Function to pack data to be sent
       
       \param[in] _fields data fields to be packed
    */
#ifdef __GXX_EXPERIMENTAL_CXX0X__
    template <typename... FIELDS>
    void pack(const FIELDS&... _fields) const {
      hd.pack(_fields... );
    }
#else
#define MACRO_IMPL(z, n, _)                                             \
    template <BOOST_PP_ENUM_PARAMS_Z(z, BOOST_PP_INC(n), typename FIELD)> \
    void pack(BOOST_PP_ENUM_BINARY_PARAMS_Z(z, BOOST_PP_INC(n), FIELD, const &_field)) { \
      hd.pack(BOOST_PP_ENUM_PARAMS_Z(z, BOOST_PP_INC(n), _field));      \
    }
  
    BOOST_PP_REPEAT(GCL_MAX_FIELDS, MACRO_IMPL, all)
#undef MACRO_IMPL
#endif

    /**
       Function to unpack received data

       \param[in] _fields data fields where to unpack data
    */
#ifdef __GXX_EXPERIMENTAL_CXX0X__
    template <typename... FIELDS>
    void unpack(const FIELDS&... _fields) const {
      hd.unpack(_fields... );
    }
#else
#define MACRO_IMPL(z, n, _)                                             \
    template <BOOST_PP_ENUM_PARAMS_Z(z, BOOST_PP_INC(n), typename FIELD)> \
    void unpack(BOOST_PP_ENUM_BINARY_PARAMS_Z(z, BOOST_PP_INC(n), FIELD, const &_field)) { \
      hd.unpack(BOOST_PP_ENUM_PARAMS_Z(z, BOOST_PP_INC(n), _field));    \
    }
    
    BOOST_PP_REPEAT(GCL_MAX_FIELDS, MACRO_IMPL, all)
#undef MACRO_IMPL
#endif

    /**
       Function to unpack received data

       \tparam array_of_fotf this should be an array of field_on_the_fly
       \param[in] fields vector with fields on the fly
    */
    template <typename T1, typename T2, template <typename> class T3>
    void pack(std::vector<field_on_the_fly<T1,T2,T3> > const& fields) {
      hd.pack(fields);
    }

    /**
       Function to unpack received data

       \tparam array_of_fotf this should be an array of field_on_the_fly
       \param[in] fields vector with fields on the fly
    */
    template <typename T1, typename T2, template <typename>class T3>
    void unpack(std::vector<field_on_the_fly<T1,T2,T3> > const& fields) {
      hd.unpack(fields);
    }

    /**
       function to trigger data exchange

       Note: when the start_exchange() + wait() combination is used, the exchange() method should not be used, and vice versa.
    */
    void exchange() {
      hd.exchange();
    }

      void post_receives() {
          hd.post_receives();
      }

      void do_sends() {
          hd.do_sends();
      }

    /**
       function to trigger data exchange initiation when using split-phase communication.

       Note: when the start_exchange() + wait() combination is used, the exchange() method should not be used, and vice versa.
    */
    void start_exchange() {
      hd.start_exchange();
    }

    /**
       function to trigger data exchange

       Note: when the start_exchange() + wait() combination is used, the exchange() method should not be used, and vice versa.
    */
    void wait() {
      hd.wait();
    }

  };
 

  template <typename layout2proc_map,
            int DIMS,
            typename Gcl_Arch = gcl_cpu,
            int version = version_manual>
  class halo_exchange_generic: public halo_exchange_generic_base<layout2proc_map,DIMS,Gcl_Arch,version> {

    typedef halo_exchange_generic_base<layout2proc_map,DIMS,Gcl_Arch,version> base_type;
    //    typedef typename layout_transform<layout_map, layout2proc_map_abs>::type layout2proc_map;

  public:
    typedef typename base_type::grid_type grid_type;

    typedef typename base_type::pattern_type pattern_type;

    template <typename DT>
    struct traits {
      static const int I=DIMS;
      typedef empty_field<DT,I> base_field;
    };

    explicit halo_exchange_generic(typename grid_type::period_type const &c, MPI_Comm comm)
      : base_type(c, comm)
    {}

    explicit halo_exchange_generic(grid_type const& g)
      : base_type(g)
    {}

  };
 

  // different traits are needed
  template <typename layout2proc_map,
            int DIMS>
  class halo_exchange_generic<layout2proc_map,DIMS,gcl_cpu,version_manual> 
    : public halo_exchange_generic_base<layout2proc_map,DIMS,gcl_cpu,version_manual> {

    //    typedef typename layout_transform<layout_map, layout2proc_map_abs>::type layout2proc_map;
    static const int version=version_manual;
    typedef gcl_cpu Gcl_Arch;

    typedef halo_exchange_generic_base<layout2proc_map,DIMS,gcl_cpu,version_manual> base_type;
  public:
    typedef typename base_type::grid_type grid_type;

    typedef typename base_type::pattern_type pattern_type;

    template <typename DT>
    struct traits {
      static const int I=DIMS;
      typedef empty_field_no_dt<I> base_field;
    };

    explicit halo_exchange_generic(typename grid_type::period_type const &c, MPI_Comm comm)
      : base_type(c, comm)
    {}

    explicit halo_exchange_generic(grid_type const& g)
      : base_type(g)
    {}

  };


  // different traits are needed
  template <typename layout2proc_map>
  class halo_exchange_generic<layout2proc_map,3,gcl_gpu,version_manual>
    : public halo_exchange_generic_base<layout2proc_map,3,gcl_gpu,version_manual> {
    static const int DIMS=3;

    //    typedef typename layout_transform<layout_map, layout2proc_map_abs>::type layout2proc_map;
    static const int version=version_manual;
    typedef gcl_gpu Gcl_Arch;
    typedef halo_exchange_generic_base<layout2proc_map,DIMS,gcl_gpu,version_manual> base_type;
  public:
    typedef typename base_type::grid_type grid_type;

    typedef typename base_type::pattern_type pattern_type;

    template <typename DT>
    struct traits {
      static const int I=DIMS;
      typedef empty_field_no_dt_gpu<I> base_field;
    };

    explicit halo_exchange_generic(typename grid_type::period_type const &c, MPI_Comm comm)
      : base_type(c, comm)
    {}

    explicit halo_exchange_generic(grid_type const& g)
      : base_type(g)
    {}

  };
}

#endif
