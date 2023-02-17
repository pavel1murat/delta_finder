///////////////////////////////////////////////////////////////////////////////
// MC datasets for the vertical beam scan
///////////////////////////////////////////////////////////////////////////////
#ifndef __init_su2020_datasets__
#define __init_su2020_datasets__

#include "Stntuple/val/stn_dataset.hh"
#include "Stntuple/val/stn_book.hh"
//-----------------------------------------------------------------------------
namespace delta_finder {
  
  void init_histogram_datasets(stn_book* Book) {
    stn_dataset_t* ds;
    hist_file_t*   hf;
//-----------------------------------------------------------------------------
//  the best way to manage multiple datasets is to arrange them alphabetically
//-----------------------------------------------------------------------------
    if (strcmp(Book->GetName(),"delta_finder") == 0) {
//-----------------------------------------------------------------------------
// delta_finder datasets
//-----------------------------------------------------------------------------
      ds = Book->NewDataset ("cele0s41b2" ,"",-1,   1000000000); 

      hf = Book->NewHistFile(ds->id(),"deltaFinder_test0","" );
      hf = Book->NewHistFile(ds->id(),"deltaFinder_test1","" );
      hf = Book->NewHistFile(ds->id(),"flagBkgHits_test0","" );
      hf = Book->NewHistFile(ds->id(),"flagBkgHits_test1","" );
      hf = Book->NewHistFile(ds->id(),"flagBkgHits_test2","" );
      hf = Book->NewHistFile(ds->id(),"flagBkgHits_test3","" );

      ds = Book->NewDataset ("pbar1s41b0" ,"",-1,   1000000000); 

      hf = Book->NewHistFile(ds->id(),"deltaFinder_test0","" );
      hf = Book->NewHistFile(ds->id(),"deltaFinder_test1","" );
      hf = Book->NewHistFile(ds->id(),"flagBkgHits_test0","" );
      hf = Book->NewHistFile(ds->id(),"flagBkgHits_test2","" );
//-----------------------------------------------------------------------------
// new version, new naming conventions, more histograms...
//-----------------------------------------------------------------------------
      ds = Book->NewDataset ("cele0b2s41r00" ,"",-1,   1000000000); 
      hf = Book->NewHistFile(ds->id(),"deltaFinder_test1","" );
//-----------------------------------------------------------------------------
// new version, new naming conventions, more histograms...
//-----------------------------------------------------------------------------
      ds = Book->NewDataset ("cele0b2s41r03" ,"",-1,   1000000000); 
      hf = Book->NewHistFile(ds->id(),"deltaFinder_test1","" );
      hf = Book->NewHistFile(ds->id(),"deltaFinder_diag" ,"" );
      hf = Book->NewHistFile(ds->id(),"flagBkgHist_test3","" );
//-----------------------------------------------------------------------------
// new version, new naming conventions, more histograms...
//-----------------------------------------------------------------------------
      ds = Book->NewDataset ("mnbs0b1s41r03" ,"",-1,   1000000000); 
      hf = Book->NewHistFile(ds->id(),"deltaFinder_test1","" );
//-----------------------------------------------------------------------------
// new version, new naming conventions, more histograms...
//-----------------------------------------------------------------------------
      ds = Book->NewDataset ("pbar1b0s41r03" ,"",-1,   1000000000); 
      hf = Book->NewHistFile(ds->id(),"deltaFinder_test1","" );
//-----------------------------------------------------------------------------
// new version, new naming conventions, more histograms...
//-----------------------------------------------------------------------------
      ds = Book->NewDataset ("cele0b2s41r04" ,"",-1,   1000000000); 
      hf = Book->NewHistFile(ds->id(),"deltaFinder_test1","" );
//-----------------------------------------------------------------------------
// new version, new naming conventions, more histograms...
//-----------------------------------------------------------------------------
      ds = Book->NewDataset ("mnbs0b1s41r04" ,"",-1,   1000000000); 
      hf = Book->NewHistFile(ds->id(),"deltaFinder_test1","" );
//-----------------------------------------------------------------------------
// new version, new naming conventions, more histograms...
//-----------------------------------------------------------------------------
      ds = Book->NewDataset ("pbar1b0s41r04" ,"",-1,   1000000000); 
      hf = Book->NewHistFile(ds->id(),"deltaFinder_test1","" );
    }
  }
}
#endif
