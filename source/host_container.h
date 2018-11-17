#pragma once
#include <memory>
#include <boost/multi_index_container.hpp>
#include <boost/multi_index/ordered_index.hpp>
#include <boost/multi_index/hashed_index.hpp>
#include <boost/multi_index/random_access_index.hpp>
#include <boost/multi_index/mem_fun.hpp>
#include "host.h"

namespace bmi = ::boost::multi_index;

struct id_tag{};
struct trial_tag{};
struct index_tag{};

typedef bmi::multi_index_container< 
  std::shared_ptr<Host>,
  bmi::indexed_by<
    bmi::hashed_unique<
      bmi::tag<id_tag>,
      bmi::const_mem_fun<Host, int, &Host::get_id>
    >,
//    bmi::hashed_non_unique<
//      bmi::tag<trial_tag>,
//      bmi::const_mem_fun<Host, std::string, &Host::get_trial_arm_name>
//    >,
    bmi::random_access<
      bmi::tag<index_tag> 
    >
  >
> HostMultiIndexContainer;