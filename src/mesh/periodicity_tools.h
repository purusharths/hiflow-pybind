// Copyright (C) 2011-2021 Vincent Heuveline
//
// HiFlow3 is free software: you can redistribute it and/or modify it under the
// terms of the European Union Public Licence (EUPL) v1.2 as published by the
// European Union or (at your option) any later version.
//
// HiFlow3 is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
// A PARTICULAR PURPOSE. See the European Union Public Licence (EUPL) v1.2 for
// more details.
//
// You should have received a copy of the European Union Public Licence (EUPL)
// v1.2 along with HiFlow3.  If not, see
// <https://joinup.ec.europa.eu/page/eupl-text-11-12>.

#ifndef PERIODICITY_TOOLS_H
#define PERIODICITY_TOOLS_H

/// \file periodicity_tools.h
/// \brief Periodicity tools.
///
/// \author Teresa Beck, Jonathan Schwegler

#include "mesh/types.h"
#include <mpi.h>

/// \brief Class MasterSlave contains relevant information for periodicity.
/// \param master Master bound
/// \param slave Slave bound
/// \param h cell width in vicinity of periodic bound
/// \param index direction of periodicity: 0 = periodicity in x-direction,
/// 1 = periodicity in y-direction, 2 = periodicity in z-direction.

namespace hiflow {
namespace mesh {

enum PeriodicBoundaryType {
  NO_PERIODIC_BDY = 0,
  MASTER_PERIODIC_BDY = 1,
  SLAVE_PERIODIC_BDY = -1
};

class MasterSlave;

template < class DataType >
std::vector< std::vector< PeriodicBoundaryType > >
get_periodicity_type(const std::vector< DataType > &coordinates,
                     const GDim gdim, const std::vector< MasterSlave > &period);

/// \brief Map points from slave boundary to master boundary. See also
/// MasterSlave
template < class DataType >
std::vector< DataType > periodify_slave_to_master(
    const std::vector< DataType > &coordinates, const GDim gdim,
    const std::vector< MasterSlave > &period, int which_period = -1);

/// \brief Map points from master boundary to slave boundary. See also
/// MasterSlave
template < class DataType >
std::vector< DataType > periodify_master_to_slave(
    const std::vector< DataType > &coordinates, const GDim gdim,
    const std::vector< MasterSlave > &period, int which_period = -1);

/// \brief Reverse the process of periodify. This function will only work
/// correctly if h of MasterSlave is chosen "good". Also the mesh should be
/// refined enough. ATTENTION: The coordinates have to form an entity.
/// To ensure this only use coordinates obtained via
/// entity.get_coordinates(std::vector<T>& coords).
template < class DataType >
std::vector< DataType > unperiodify(const std::vector< DataType > &coordinates,
                                    const GDim gdim,
                                    const std::vector< MasterSlave > &period);

template < class DataType >
std::vector< std::vector< DataType > >
unperiodify(const std::vector< std::vector< DataType > > &coordinates,
            const GDim gdim, const std::vector< MasterSlave > &period);

} // namespace mesh
} // namespace hiflow
#endif
