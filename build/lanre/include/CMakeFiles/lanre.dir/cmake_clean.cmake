file(REMOVE_RECURSE
  "../../../lib/liblanre.dylib"
  "../../../lib/liblanre.pdb"
)

# Per-language clean rules from dependency scanning.
foreach(lang )
  include(CMakeFiles/lanre.dir/cmake_clean_${lang}.cmake OPTIONAL)
endforeach()
