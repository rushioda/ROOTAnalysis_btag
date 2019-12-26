################################################################################
# Build configuration for MyAnalysisExampleROOT
################################################################################

# Declare the name of the package:
atlas_subdir( MyAnalysisExampleROOT )

atlas_depends_on_subdirs( PUBLIC
			  PhysicsAnalysis/D3PDTools/EventLoop
			  Event/xAOD/xAODBase
			  Control/xAODRootAccess
			  Event/xAOD/xAODEventInfo
        Event/xAOD/xAODMuon
        Event/xAOD/xAODTruth
	  		Control/AthToolSupport/AsgTools
        PhysicsAnalysis/MCTruthClassifier
        InnerDetector/InDetRecTools/InDetTrackSelectionTool
)

# External dependencies:
find_package( ROOT COMPONENTS Core RIO Hist Physics Boost )

# Build a dictionary for the library
atlas_add_root_dictionary( MyAnalysisExampleROOT _dictionarySource
			   ROOT_HEADERS Root/LinkDef.h
			   EXTERNAL_PACKAGES ROOT )

# Build a shared library:
atlas_add_library( MyAnalysisExampleROOT
   MyAnalysisExampleROOT/*.h Root/*.h Root/*.cxx ${_dictionarySource}
   INCLUDE_DIRS ${ROOT_INCLUDE_DIRS}
   PUBLIC_HEADERS MyAnalysisExampleROOT
   LINK_LIBRARIES ${ROOT_LIBRARIES} AsgTools EventLoop xAODBase xAODRootAccess xAODEventInfo xAODMuon xAODTruth MCTruthClassifierLib EventLoopAlgs InDetTrackSelectionToolLib xAODJet)

# Executable
atlas_add_executable( testRun
		      util/testRun.cxx
		      INCLUDE_DIRS ${ROOT_INCLUDE_DIRS}
		      LINK_LIBRARIES ${ROOT_LIBRARIES} AsgTools EventLoop xAODBase xAODRootAccess xAODEventInfo xAODMuon xAODTruth MyAnalysisExampleROOT )
