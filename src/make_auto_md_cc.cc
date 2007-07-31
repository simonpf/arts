/* Copyright (C) 2000-2007 Stefan Buehler <sbuehler@ltu.se>

   This program is free software; you can redistribute it and/or modify it
   under the terms of the GNU General Public License as published by the
   Free Software Foundation; either version 2, or (at your option) any
   later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307,
   USA. */

#include "arts.h"
#include "token.h"
#include "array.h"
#include "file.h"
#include "auto_wsv.h"
#include "methods.h"
#include "wsv_aux.h"
#include "agenda_record.h"

/* Adds commas and indentation to parameter lists. */
void align(ofstream& ofs, bool& is_first_parameter, const String& indent)
{
  // Add comma and line break, if not first element:
  if (is_first_parameter)
    is_first_parameter = false;
  else
    {
      ofs << ",\n";
      // Make proper indentation:
      ofs << indent;
    }
}

int main()
{
  try
    {
      // Make the global data visible:
      extern Array<MdRecord> md_data;
      extern const ArrayOfString wsv_group_names;
      extern const Array<WsvRecord> wsv_data;

      // Initialize method data.
      define_md_data_raw();

      // Initialize the wsv group name array:
      define_wsv_group_names();

      // Expand supergeneric methods:
      expand_md_data_raw_to_md_data();

      // Initialize wsv data.
      define_wsv_data();
  

      const Index n_md  = md_data.nelem();
      const Index n_wsv = wsv_data.nelem();

      // For safety, check if n_wsv and N_WSV have the same value. If not, 
      // then the file wsv.h is not up to date.
      if (N_WSV != n_wsv)
        {
          cout << "The file wsv.h is not up to date!\n";
          cout << "(N_WSV = " << N_WSV << ", n_wsv = " << n_wsv << ")\n";
          cout << "Make wsv.h first. Check if Makefile is correct.\n";
          return 1;
        }

      // Write auto_md.cc:
      // -----------
      ofstream ofs;
      open_output_file(ofs,"auto_md.cc");
  
      ofs << "// This file was generated automatically by make_auto_md_cc.cc.\n";
      ofs << "// DO NOT EDIT !\n";
      ofs << "// Generated: "
          << __DATE__ << ", "
          << __TIME__ << "\n\n";

      ofs << "#include \"arts.h\"\n"
          << "#include \"make_array.h\"\n"
          << "#include \"auto_md.h\"\n"
          << "#include \"auto_wsv_groups.h\"\n"
          << "#include \"wsv_aux.h\"\n"
          << "#include \"m_copy.h\"\n"
          << "#include \"m_general.h\"\n"
          << "#include \"m_ignore.h\"\n"
          << "#include \"m_xml.h\"\n"
          << "#include \"m_basic_types.h\"\n"
          << "#include \"agenda_record.h\"\n"
          << "\n";

      // Declare wsv_data:
      ofs << "// Other wsv data:\n"
          << "extern const Array<WsvRecord> wsv_data;\n\n";


      // Write all get-away functions:
      // -----------------------------
      for (Index i=0; i<n_md; ++i)
        {
          const MdRecord& mdd = md_data[i];

          // This is needed to flag the first function parameter, which 
          // needs no line break before being written:
          bool is_first_parameter = true;
          // The String indent is needed to achieve the correct
          // indentation of the functin parameters:
          String indent = String(mdd.Name().nelem()+3,' ');;
          
          // There are four lists of parameters that we have to
          // write. 
          ArrayOfIndex  vo=mdd.Output();   // Output 
          ArrayOfIndex  vi=mdd.Input();    // Input
          ArrayOfIndex  vgo=mdd.GOutput();   // Generic Output 
          ArrayOfIndex  vgi=mdd.GInput();    // Generic Input
          // vo and vi contain handles of workspace variables, 
          // vgo and vgi handles of workspace variable groups.

          // Check, if some workspace variables are in both the
          // input and the output list, and erase those from the input 
          // list:
          for (ArrayOfIndex::const_iterator j=vo.begin(); j<vo.end(); ++j)
            for (ArrayOfIndex::iterator k=vi.begin(); k<vi.end(); ++k)
              if ( *j == *k )
                {
                  //              erase_vector_element(vi,k);
                  k = vi.erase(k) - 1;
                  // We need the -1 here, otherwise due to the
                  // following increment we would miss the element
                  // behind the erased one, which is now at the
                  // position of the erased one.
                }

          // There used to be a similar block here for the generic
          // input/output variables. However, this was a mistake. For
          // example, if a method has a vector as generic input and a
          // vector as generic output, this does not mean that it is
          // the same vector!

            {

              String ws, mr;

              // Use parameter name only if it is used inside the function
              // to avoid warnings
              ws = " ws";
              if (!vo.nelem () && !vi.nelem () && !vgo.nelem () && !vgi.nelem ())
              {
                ws = "";
              }

              // Use parameter name only if it is used inside the function
              // to avoid warnings
              if ( vgo.nelem () || vgi.nelem ()
                   || mdd.Keywords().nelem() || mdd.AgendaMethod())
                {
                  mr = " mr";
                }

              if ( mdd.Supergeneric() )
                {
                  ofs << "void " << mdd.Name()
                    << "_sg_" << wsv_group_names[mdd.ActualGroup()]
                    << "_g(Workspace&" << ws
                    << ", const MRecord&" << mr << ")\n"
                    << "{\n";
                }
              else
                {
                  ofs << "void " << mdd.Name()
                    << "_g(Workspace&" << ws
                    << ", const MRecord&" << mr << ")\n"
                    << "{\n";
                }
            }

          // Define generic output pointers
          for (Index j=0; j<vgo.nelem(); ++j)
            {
              // Check by assert whether the group identifier is too
              // large to correspond to a group. This can easily
              // happen if somebody puts a variable identifier instead
              // of a group identifier in the argument of GOUTPUT:
              assert( vgo[j] < wsv_group_names.nelem() );

              ofs << "  " << wsv_group_names[vgo[j]]
                  << " *GO" << j << " = (" << wsv_group_names[vgo[j]] << " *)ws[mr.Output()[" << j
                  << "]];\n";
            }

          // Define generic input pointers
          for (Index j=0; j<vgi.nelem(); ++j)
            {
              // Check by assert whether the group identifier is too
              // large to correspond to a group. This can easily
              // happen if somebody puts a variable identifier instead
              // of a group identifier in the argument of GINPUT:
              assert( vgi[j] < wsv_group_names.nelem() );

              ofs << "  " << wsv_group_names[vgi[j]]
                  << " *GI" << j << " = (" << wsv_group_names[vgi[j]] << " *)ws[mr.Input()[" << j
                  << "]];\n";
            }

          ofs << "  " << mdd.Name() << "(";

          // Write the Output workspace variables:
          for (Index j=0; j<vo.nelem(); ++j)
            {
              // Add comma and line break, if not first element:
              align(ofs,is_first_parameter,indent);

              //ofs << "ws." << wsv_data[vo[j]].Name();
              ofs << "*((" << wsv_group_names[wsv_data[vo[j]].Group()] <<  " *)ws["
                << vo[j] << "])";
            }

          // Write the Generic output workspace variables:
          for (Index j=0; j<vgo.nelem(); ++j)
            {
              // Add comma and line break, if not first element:
              align(ofs,is_first_parameter,indent);

              ofs << "*GO" << j;
            }

          // Write the Generic output workspace variable names:
          for (Index j=0; j<vgo.nelem(); ++j)
            {
              // Add comma and line break, if not first element:
              align(ofs,is_first_parameter,indent);

              ofs << "wsv_data[mr.Output()["
                  << j
                  << "]].Name()";
            }

          // Write the Input workspace variables:
          for (Index j=0; j<vi.nelem(); ++j)
            {
              // Add comma and line break, if not first element:
              align(ofs,is_first_parameter,indent);

              //ofs << "ws." << wsv_data[vi[j]].Name();
              ofs << "*((" << wsv_group_names[wsv_data[vi[j]].Group()] <<  " *)ws["
                << vi[j] << "])";
            }

          // Write the Generic input workspace variables:
          for (Index j=0; j<vgi.nelem(); ++j)
            {
              // Add comma and line break, if not first element:
              align(ofs,is_first_parameter,indent);

              ofs << "*GI" << j;
            }

          // Write the Generic input workspace variable names:
          for (Index j=0; j<vgi.nelem(); ++j)
            {
              // Add comma and line break, if not first element:
              align(ofs,is_first_parameter,indent);

              ofs << "wsv_data[mr.Input()["
                  << j
                  << "]].Name()";
            }

          // Write the control parameters:
          {
            // The mr parameters look all the same (mr[i]), so we just
            // need to know the number of them: 
            Index n_mr = mdd.Keywords().nelem();

            for (Index j=0; j!=n_mr; ++j)
              {
                // Add comma and line break, if not first element:
                align(ofs,is_first_parameter,indent);

                ofs << "mr.Values()[" << j << "]";
              }
          }

          // Write the agenda, if there is one.
          if ( mdd.AgendaMethod() )
            {
              align(ofs,is_first_parameter,indent);
              ofs << "mr.Tasks()";
            }

          ofs << ");\n";
          ofs << "}\n\n";
        }

      // Add getaways, the array that hold pointers to the getaway functions:
      {
        String indent = "     ";
        bool is_first_parameter = true;

        ofs << "// The array holding the pointers to the getaway functions.\n"
            << "void (*getaways[])(Workspace&, const MRecord&)\n"
            << "  = {";
        for (Index i=0; i<n_md; ++i)
          {
            const MdRecord& mdd = md_data[i];
          
            // Add comma and line break, if not first element:
            align(ofs,is_first_parameter,indent);

            if ( mdd.Supergeneric() )
              {
                ofs << mdd.Name()
                    << "_sg_" << wsv_group_names[mdd.ActualGroup()]
                    << "_g";
              }
            else
              {
                ofs << mdd.Name() << "_g";
              }
          }
        ofs << "};\n\n";
      }

      // Create implementation of the agenda wrappers

      // Initialize agenda data.
      define_wsv_map ();
      define_agenda_data ();

      extern const Array<AgRecord> agenda_data;
      for (Index i = 0; i < agenda_data.nelem (); i++)
        {
          const AgRecord& agr = agenda_data[i];
          const ArrayOfIndex& ago = agr.Output ();
          const ArrayOfIndex& agi = agr.Input ();
          ostringstream ain_push_os, ain_pop_os;
          ostringstream aout_push_os, aout_pop_os;

          write_agenda_wrapper_header (ofs, agr);

          ofs << "\n";
          ofs << "{\n";
          ofs << "  extern Workspace workspace;\n";
          if (ago.nelem () || agi.nelem ())
            {
              ofs << "  extern map<String, Index> AgendaMap;\n"
                << "  extern const Array<AgRecord> agenda_data;\n"
                << "\n"
                << "  const AgRecord& agr =\n"
                << "    agenda_data[AgendaMap.find (input_agenda.name ())->second];\n"
                << "\n";
            }
          if (ago.nelem ())
            {
              for (Index j = 0; j < ago.nelem (); j++)
                {
                  // Mark agenda output-only variables as uninitialized
                  ArrayOfIndex::const_iterator it = agi.begin ();
                  while (it != agi.end () && *it != ago[j]) it++;
                  if (it == agi.end ())
                    {
                      aout_push_os << "  workspace.push_uninitialized (aout[" << j << "], "
                        << "(void *)&" << wsv_data[ago[j]].Name () << ");\n";
                    }
                  else
                    {
                      aout_push_os << "  workspace.push (aout[" << j << "], "
                        << "(void *)&" << wsv_data[ago[j]].Name () << ");\n";
                    }
                  aout_pop_os << "  workspace.pop (aout[" << j << "]);\n";
                }
            }
          if (agi.nelem ())
            {
              for (Index j = 0; j < agi.nelem (); j++)
                {
                  // Ignore Input parameters that are also output
                  ArrayOfIndex::const_iterator it = ago.begin ();
                  while (it != ago.end () && *it != agi[j]) it++;
                  if (it == ago.end ())
                    {
                      ain_push_os << "  workspace.push (ain[" << j << "], "
                        << "(void *)&" << wsv_data[agi[j]].Name () << ");\n";
                      ain_pop_os << "  workspace.pop (ain[" << j << "]);\n";
                    }
                }
            }

          if (aout_push_os.str().length())
            {
              ofs << "  const ArrayOfIndex& aout = agr.Output ();\n";
              ofs << aout_push_os.str () << "\n";
            }
          if (ain_push_os.str().length())
            {
              ofs << "  const ArrayOfIndex& ain = agr.Input ();\n";
              ofs << ain_push_os.str () << "\n";
            }

          ofs << "  const ArrayOfIndex& outputs_to_push = input_agenda.get_output2push();\n"
              << "  const ArrayOfIndex& outputs_to_dup = input_agenda.get_output2dup();\n"
              << "\n"
              << "  for (ArrayOfIndex::const_iterator it = outputs_to_push.begin ();\n"
              << "       it != outputs_to_push.end (); it++)\n"
              << "  { workspace.push (*it, NULL); }\n"
              << "\n"
              << "  for (ArrayOfIndex::const_iterator it = outputs_to_dup.begin ();\n"
              << "       it != outputs_to_dup.end (); it++)\n"
              << "  { workspace.duplicate (*it); }\n"
              << "\n";

          ofs << "  input_agenda.execute (silent);\n\n";

          ofs << "  for (ArrayOfIndex::const_iterator it = outputs_to_push.begin ();\n"
              << "       it != outputs_to_push.end (); it++)\n"
              << "    { workspace.pop_free (*it); }\n"
              << "\n"
              << "  for (ArrayOfIndex::const_iterator it = outputs_to_dup.begin ();\n"
              << "       it != outputs_to_dup.end (); it++)\n"
              << "    { workspace.pop_free (*it); }\n\n";

          if (aout_pop_os.str().length())
            {
              ofs << aout_pop_os.str () << "\n";
            }
          if (ain_pop_os.str().length())
            {
              ofs << ain_pop_os.str () << "\n";
            }

          ofs << "}\n\n";
        }
    }
  catch (exception x)
    {
      cout << "Something went wrong. Message text:\n";
      cout << x.what() << '\n';
      return 1;
    }

  return 0;
}
