//
// Created by daniel on 04/10/22.
//

#ifndef ABMCMC_BOOSTSERIALIZATION_H
#define ABMCMC_BOOSTSERIALIZATION_H

namespace boost {
    namespace serialization {

        template<class Archive, class T>
        inline void serialize(
                Archive & ar,
                std::complex< T > & t,
                const unsigned int file_version
        ){
            boost::serialization::split_free(ar, t, file_version);
        }

        template<class Archive, class T>
        inline void save(
                Archive & ar,
                std::complex< T > const & t,
                const unsigned int /* file_version */
        ){
            const T re = t.real();
            const T im = t.imag();
            ar << boost::serialization::make_nvp("real", re);
            ar << boost::serialization::make_nvp("imag", im);
        }

        template<class Archive, class T>
        inline void load(
                Archive & ar,
                std::complex< T >& t,
                const unsigned int /* file_version */
        ){
            T re;
            T im;
            ar >> boost::serialization::make_nvp("real", re);
            ar >> boost::serialization::make_nvp("imag", im);
            t = std::complex< T >(re,im);
        }

// specialization of serialization traits for complex
        template <class T>
        struct is_bitwise_serializable<std::complex< T > >
        : public is_bitwise_serializable< T > {};

    template <class T>
    struct implementation_level<std::complex< T > >
    : mpl::int_<object_serializable> {} ;

// treat complex just like builtin arithmetic types for tracking
template <class T>
struct tracking_level<std::complex< T > >
: mpl::int_<track_never> {} ;

} // serialization
} // namespace boost


#endif //ABMCMC_BOOSTSERIALIZATION_H
