/*!
  \file   groups.cc
  \brief  Defines workspace variable groups.

  If you want to add new workspace variable groups you have to do it
  in this file. This is used by the program make_wsv_group_h to
  generate the header file wsv_group.h

  \author Stefan Buehler
  \date   2000-08-04 */

#include "groups.h"

#include <map>

#include "array.h"
#include "arts.h"
#include "mystring.h"

/*! The names associated with Wsv groups as Strings.
  See function define_wsv_groups for more information. */
namespace global_data {
ArrayOfGroupRecord wsv_groups;
map<String, Index> WsvGroupMap;
}  // namespace global_data

void define_wsv_group_map() {
  using global_data::wsv_groups;
  using global_data::WsvGroupMap;
  for (Index i = 0; i < wsv_groups.nelem(); ++i) {
    WsvGroupMap[wsv_groups[i]] = i;
  }
}

//! Define the array of workspace variable group names.
/*!
  This defines the global variable wsv_groups. It is used in two
  different programs:

  1. In arts.

  2. In make_wsv_group_h.

  \author Stefan Buehler
  \date   2000-08-04
*/
void define_wsv_groups() {
  using global_data::wsv_groups;

  //--------------------< Build the group names array >--------------------
  // Initialize to empty, just in case.
  wsv_groups.resize(0);

  wsv_groups.emplace_back(
      "AbsorptionLines",
      "Contains line-by-line absorption information for a number of related absorption lines");

  wsv_groups.emplace_back(
      "Agenda", "Describes a set of function calls and variable definitions");

  wsv_groups.emplace_back(
      "Any",
      "Meta type for when methods can take any argument (avoid manual use)");

  wsv_groups.emplace_back("ArrayOfAbsorptionLines",
                          "A list of *AbsorptionLines*");

  wsv_groups.emplace_back("ArrayOfArrayOfAbsorptionLines",
                          "A list of *ArrayOfAbsorptionLines*");

  wsv_groups.emplace_back("ArrayOfAgenda", "A list of *Agenda*");

  wsv_groups.emplace_back("ArrayOfArrayOfGriddedField1",
                          "A list of *ArrayOfGriddedField1*");

  wsv_groups.emplace_back("ArrayOfArrayOfGriddedField2",
                          "A list of *ArrayOfGriddedField2*");

  wsv_groups.emplace_back("ArrayOfArrayOfGriddedField3",
                          "A list of *ArrayOfGriddedField3*");

  wsv_groups.emplace_back("ArrayOfArrayOfIndex", "A list of *ArrayOfIndex*");

  wsv_groups.emplace_back("ArrayOfArrayOfMatrix", "A list of *ArrayOfMatrix*");

  wsv_groups.emplace_back("ArrayOfPpath", "A list of *Ppath*");

  wsv_groups.emplace_back("ArrayOfArrayOfScatteringMetaData",
                          "A list of *ArrayOfScatteringMetaData*");

  wsv_groups.emplace_back("ArrayOfArrayOfSingleScatteringData",
                          "A list of *ArrayOfSingleScatteringData*");

  wsv_groups.emplace_back("ArrayOfArrayOfSpeciesTag",
                          "A list of *ArrayOfSpeciesTag*");

  wsv_groups.emplace_back("ArrayOfArrayOfString", "A list of *ArrayOfString*");

  wsv_groups.emplace_back("ArrayOfArrayOfTensor3",
                          "A list of *ArrayOfTensor3*");

  wsv_groups.emplace_back("ArrayOfArrayOfTensor6",
                          "A list of *ArrayOfTensor6*");

  wsv_groups.emplace_back("ArrayOfArrayOfTime", "A list of *ArrayOfTime*");

  wsv_groups.emplace_back("ArrayOfArrayOfVector", "A list of *ArrayOfVector*");

  wsv_groups.emplace_back("ArrayOfAtmPoint", "A list of *AtmPoint*");

  wsv_groups.emplace_back("ArrayOfCIARecord", "A list of *CIARecord*");

  wsv_groups.emplace_back("ArrayOfGriddedField1", "A list of *GriddedField1*");

  wsv_groups.emplace_back("ArrayOfGriddedField2", "A list of *GriddedField2*");

  wsv_groups.emplace_back("ArrayOfGriddedField3", "A list of *GriddedField3*");

  wsv_groups.emplace_back("ArrayOfGriddedField4", "A list of *GriddedField4*");

  wsv_groups.emplace_back("ArrayOfIndex", "A list of *Index*");

  wsv_groups.emplace_back("ArrayOfJacobianTarget",
                          "A list of *JacobianTarget*");

  wsv_groups.emplace_back("ArrayOfMatrix", "A list of *Matrix*");

  wsv_groups.emplace_back("ArrayOfQuantumIdentifier",
                          "A list of *QuantumIdentifier*");

  wsv_groups.emplace_back("ArrayOfRetrievalQuantity",
                          "A list of retrieval quantitities");

  wsv_groups.emplace_back("ArrayOfScatteringMetaData",
                          "A list of *ScatteringMetaData*");

  wsv_groups.emplace_back("ArrayOfSingleScatteringData",
                          "A list of *SingleScatteringData*");

  wsv_groups.emplace_back("ArrayOfSpeciesTag", R"--(A list of species tags

These tags include the species and a lot of optional information
about the isotopologue, the absorption scheme, and the frequency limits)--");

  wsv_groups.emplace_back("ArrayOfSparse", "A list of *Sparse*");

  wsv_groups.emplace_back("ArrayOfSun", "A list of sun");

  wsv_groups.emplace_back("ArrayOfString", "A list of *String*");

  wsv_groups.emplace_back("ArrayOfTelsemAtlas", "A list of *TelsemAtlas*");

  wsv_groups.emplace_back("ArrayOfTensor3", "A list of *Tensor3*");

  wsv_groups.emplace_back("ArrayOfTensor4", "A list of *Tensor4*");

  wsv_groups.emplace_back("ArrayOfTensor5", "A list of *Tensor5*");

  wsv_groups.emplace_back("ArrayOfTensor6", "A list of *Tensor6*");

  wsv_groups.emplace_back("ArrayOfTensor7", "A list of *Tensor7*");

  wsv_groups.emplace_back("ArrayOfTime", "A list of *Time*");

  wsv_groups.emplace_back("ArrayOfVector", "A list of *Vector*");

  wsv_groups.emplace_back("ArrayOfXsecRecord",
                          R"--(A list of cross-section records

These cross-section records contains information about the valid temperature and
pressure ranges as well as well as the fitting coefficients used to compute
and interpolate the cross-section to other temperatures and pressures)--");

  wsv_groups.emplace_back("AtmField", R"--(An atmospheric field)--");

  wsv_groups.emplace_back("AtmPoint", R"--(An atmospheric point)--");

  wsv_groups.emplace_back(
      "CIARecord",
      R"--(Contains information to compute collision induced absorption for a pair of species

Holds an the record data in a gridded field with grids of temperature and frequency in
units of m^5 molec^(-2) )--");

  wsv_groups.emplace_back("CallbackFunction",
                          "Used to inject custom code into *Agenda*");

  wsv_groups.emplace_back("CovarianceMatrix", "Contains the covariance matrix");

  wsv_groups.emplace_back("GasAbsLookup", R"--(An absorption lookup table

This class holds an absorption lookup table, as well as all
information that is necessary to use the table to extract
absorption)--");

  wsv_groups.emplace_back("GridPos", "A position in a grid");

  wsv_groups.emplace_back("GriddedField1",
                          R"--(A 1 dimensional gridded set of *Numeric* data

The grid is 1 *Vector* or *ArrayOfString*

Both the data and the grid may be named)--");

  wsv_groups.emplace_back("GriddedField2",
                          R"--(A 2 dimensional gridded set *Numeric* data

The grid is a combination of 2 *Vector* and/or *ArrayOfString*

Both the data and the grid may be named)--");

  wsv_groups.emplace_back("GriddedField3",
                          R"--(A 3 dimensional gridded set of *Numeric* data

The grid is a combination of 3 *Vector* and/or *ArrayOfString*

Both the data and the grid may be named)--");

  wsv_groups.emplace_back("GriddedField4",
                          R"--(A 4 dimensional gridded set of *Numeric* data

The grid is a combination of 4 *Vector* and/or *ArrayOfString*

Both the data and the grid may be named)--");

  wsv_groups.emplace_back("GriddedField5",
                          R"--(A 5 dimensional gridded set  of *Numeric* data

The grid is a combination of 5 *Vector* and/or *ArrayOfString*

Both the data and the grid may be named)--");

  wsv_groups.emplace_back("GriddedField6",
                          R"--(A 6 dimensional gridded set of *Numeric* data

The grid is a combination of 6 *Vector* and/or *ArrayOfString*

Both the data and the grid may be named)--");

  wsv_groups.emplace_back("HitranRelaxationMatrixData",
                          "Wraps data required to use Hitran line mixing");

  wsv_groups.emplace_back("Index", "A 64 bit signed integer type");

  wsv_groups.emplace_back(
      "JacobianTarget", "A single target if a partial derivative computation");

  wsv_groups.emplace_back(
      "MapOfErrorCorrectedSuddenData",
      R"--(A map of data required for computing the error-corrected-sudden relaxation matrix

This map contains a list of an underlying data type.  This underlying data type contains a
*QuantumIdentifier* and a list of species dependent computational data for various components
required to compute the relaxation matrix

If there is no identifier or species avaialable, default values that approximates a diagonal
relaxation matrix are set)--");

  wsv_groups.emplace_back("MCAntenna", "An Antenna object used by *MCGeneral*");

  wsv_groups.emplace_back("Matrix", "A 2 dimensional array of *Numeric*");

  wsv_groups.emplace_back("Numeric", "IEEE 754 binary64 floating point number");

  wsv_groups.emplace_back("Ppath", "Describes a propagation path");

  wsv_groups.emplace_back("PredefinedModelData",
                          R"--(Contains any data required for a predefined model)--");

  wsv_groups.emplace_back("QuantumIdentifier",
                          R"--(An ID for an absorption species state

It contains information about the species and a set of quantum numbers
and can thus be used to identify one of the following:
1) a species
2) an isotopologue of a species
3) an absorption band of an isotopologue
4) an absorption line of an isotopologue
5) the energy level of absorption band(s) of an isotopologue
6) the energy level of absorption line(s) of an isotopologue)--");

  wsv_groups.emplace_back("Rational",
                          "Holds a rational number as two *Index* n / d");

  wsv_groups.emplace_back("ScatteringMetaData",
                          "Holds meta data about the scattering");

  wsv_groups.emplace_back("SingleScatteringData",
                          "Holds single scattering data");

  wsv_groups.emplace_back("Sparse", "A sparse version of *Matrix*");

  wsv_groups.emplace_back(
      "SpeciesIsotopologueRatios",
      "Contains a list of isotopologue ratios for all defined species");

  wsv_groups.emplace_back("String", "Basic string type");

  wsv_groups.emplace_back("SurfaceField",
                          R"--(A surface field that keeps relevant surface parameters)--");

  wsv_groups.emplace_back("SurfacePoint",
                          R"--(A surface point, to keep point values for the surface,
including the local normal vector)--");

  wsv_groups.emplace_back("TelsemAtlas", R"--(A telsem atlas

Represents a Telsem2 atlas containing land surface microwave emissivities.
Since the Atlas contains emissivities only for land surfaces, the data is
stored in a sparse format.
 
The emissivities are represented on an equal area grid and numbered
sequentially starting with the first latitude band at -90 degrees and
moving up to 90 degrees.

The correspondance array contains the data indices for each cellnumber
if it is contained in the Atlas and NAN otherwise.)--");

  wsv_groups.emplace_back("Tensor3", "A 3 dimensional array of *Numeric*");

  wsv_groups.emplace_back("Tensor4", "A 4 dimensional array of *Numeric*");

  wsv_groups.emplace_back("Tensor5", "A 5 dimensional array of *Numeric*");

  wsv_groups.emplace_back("Tensor6", "A 6 dimensional array of *Numeric*");

  wsv_groups.emplace_back("Tensor7", "A 7 dimensional array of *Numeric*");

  wsv_groups.emplace_back("Timer", "Represents a clock");

  wsv_groups.emplace_back("Time", R"(Represents a time stamp in the format:
"YEAR-MONTH-DAY HOUR:MINUTE:SECOND", e.g., "2023-03-06 14:32:35.35"

Note that most direct user input of a Time accepts a string as above to
represent the time stamp.
)");

  wsv_groups.emplace_back(
      "TessemNN", "Data required by TESSEM to calculate surface emissivity");

  wsv_groups.emplace_back("Vector", "A 1 dimensional array of *Numeric*");

  wsv_groups.emplace_back("VibrationalEnergyLevels", "A map of vibrational energy levels for NLTE calculations");
  
  // rtepack types
  wsv_groups.emplace_back("Propmat", R"--(A single propagation matrix.

Due to the properties of a propagation matrix, only 7 independents need be stored.
The propagation matrix is thus represented as:
  A  B  C  D
  B  A  U  V
  C -U  A  W
  D -V -W  A

This type is related to *Stokvec* in that its first 4 elements are the same as
the first 4 elements of *Stokvec* for pure clearsky radiative transfers.

This type is also related to *Muelmat* because it is computed often as the exponent
of this term multiplied by a negative distance.
)--");

wsv_groups.emplace_back("Muelmat", "A single Mueller 4x4 matrix.");

wsv_groups.emplace_back("Stokvec", "A single Stokes vector (of length 4).");

wsv_groups.emplace_back("PropmatVector", "A vector of *Propmat*.");

wsv_groups.emplace_back("MuelmatVector", "A vector of *Muelmat*.");

wsv_groups.emplace_back("StokvecVector", "A vector of *Stokvec*.");

wsv_groups.emplace_back("PropmatMatrix", "A matrix of *Propmat*.");

wsv_groups.emplace_back("MuelmatMatrix", "A matrix of *Muelmat*.");

wsv_groups.emplace_back("StokvecMatrix", "A matrix of *Stokvec*.");

wsv_groups.emplace_back("ArrayOfPropmatVector", "A list of *PropmatVector*.");

wsv_groups.emplace_back("ArrayOfMuelmatVector", "A list of *MuelmatVector*.");

wsv_groups.emplace_back("ArrayOfStokvecVector", "A list of *StokvecVector*.");

wsv_groups.emplace_back("ArrayOfPropmatMatrix", "A list of *PropmatMatrix*.");

wsv_groups.emplace_back("ArrayOfMuelmatMatrix", "A list of *MuelmatMatrix*.");

wsv_groups.emplace_back("ArrayOfStokvecMatrix", "A list of *StokvecMatrix*.");

wsv_groups.emplace_back("ArrayOfArrayOfPropmatVector", "A list of *ArrayOfPropmatVector*.");

wsv_groups.emplace_back("ArrayOfArrayOfMuelmatVector", "A list of *ArrayOfMuelmatVector*.");

wsv_groups.emplace_back("ArrayOfArrayOfStokvecVector", "A list of *ArrayOfStokvecVector*.");

wsv_groups.emplace_back("ArrayOfArrayOfPropmatMatrix", "A list of *ArrayOfPropmatMatrix*.");

wsv_groups.emplace_back("ArrayOfArrayOfMuelmatMatrix", "A list of *ArrayOfMuelmatMatrix*.");

wsv_groups.emplace_back("ArrayOfArrayOfStokvecMatrix", "A list of *ArrayOfStokvecMatrix*.");

  std::sort(wsv_groups.begin(), wsv_groups.end(), [](auto& a, auto& b) {
    return a.name < b.name;
  });

  define_wsv_group_map();
}

Index get_wsv_group_id(const String& name) {
  using global_data::WsvGroupMap;
  auto it = WsvGroupMap.find(name);
  if (it == WsvGroupMap.end()) return -1;
  return it->second;
}

void get_wsv_group_ids(ArrayOfIndex& ids, String name) {
  ids.resize(0);

  Index pos = 0;
  while (pos < name.nelem()) {
    switch (name[pos]) {
      case ' ':
      case '\r':
      case '\t':
      case '#':
        name.erase(pos, 1);
        break;
      default:
        pos++;
    }
  }

  pos = 0;
  Index prev = 0;
  while (pos < name.nelem()) {
    while (pos < name.nelem() && name[pos] != ',') pos++;
    Index id = get_wsv_group_id(name.substr(prev, pos - prev));
    if (id == -1) {
      ids.resize(0);
      return;
    }
    ids.push_back(id);
    pos++;
    prev = pos;
  }
}

bool is_agenda_group_id(const Index group) {
  return (group == get_wsv_group_id("Agenda") ||
          group == get_wsv_group_id("ArrayOfAgenda"));
}

String get_array_groups_as_string(bool basetype_is_group,
                                  bool return_basetype_only) {
  using global_data::wsv_groups;
  String arraygroups;

  bool first = true;
  for (Index i = 0; i < wsv_groups.nelem(); i++) {
    if (wsv_groups[i].name.substr(0, String("ArrayOf").length()) == "ArrayOf") {
      const String basetype = wsv_groups[i].name.substr(
          String("ArrayOf").length(), wsv_groups[i].name.length());
      bool basetype_exists = (get_wsv_group_id(basetype) != -1);

      if (return_basetype_only) {
        // Return only the basetype of the array,
        // skip arrays whose basetype is not a WSV group
        if (basetype_exists) {
          if (!first)
            arraygroups += ", ";
          else
            first = false;
          arraygroups += basetype;
        }
      } else {
        if (!basetype_is_group || (basetype_is_group && basetype_exists)) {
          if (!first)
            arraygroups += ", ";
          else
            first = false;
          arraygroups += wsv_groups[i];
        }
      }
    }
  }
  return arraygroups;
}
