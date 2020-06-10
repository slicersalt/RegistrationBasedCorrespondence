cmake_minimum_required(VERSION 3.13.4)

project(RegistrationBasedCorrespondence)

#-----------------------------------------------------------------------------
# Extension meta-information
set(EXTENSION_HOMEPAGE "https://www.slicer.org/wiki/Documentation/Nightly/Extensions/RegistrationBasedCorrespondence")
set(EXTENSION_CATEGORY "Shape Creation")
set(EXTENSION_CONTRIBUTORS "David Allemang (Kitware)")
set(EXTENSION_DESCRIPTION "Generate correspondence of consistent non-spherical topologies")
set(EXTENSION_ICONURL "http://www.example.com/Slicer/Extensions/RegistrationBasedCorrespondence.png")
set(EXTENSION_SCREENSHOTURLS "http://www.example.com/Slicer/Extensions/RegistrationBasedCorrespondence/Screenshots/1.png")
set(EXTENSION_DEPENDS "NA") # Specified as a list or "NA" if no dependencies

#-----------------------------------------------------------------------------
# Extension dependencies
find_package(Slicer REQUIRED)
include(${Slicer_USE_FILE})

#-----------------------------------------------------------------------------
# Extension modules
add_subdirectory(RegistrationBasedCorrespondence)
## NEXT_MODULE

#-----------------------------------------------------------------------------
include(${Slicer_EXTENSION_GENERATE_CONFIG})
include(${Slicer_EXTENSION_CPACK})
