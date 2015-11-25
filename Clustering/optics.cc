/* -------------------------------------------------------------------------------- */
/*                                       OPTICS                                     */
/* -------------------------------------------------------------------------------- */
/*   (c) University of Munich, Database Group                                       */
/*   written by M. Breunig, breunig@dbs.informatik.uni-muenchen.de                  */
/* -------------------------------------------------------------------------------- */
/* This program implements a very simple version of the OPTICS algorithm published  */
/* in "Ankerst M., Breunig M. M., Kriegel H.-P., Sander J.: OPTICS: Ordering Points */
/* To Identify the Clustering Structure, Proc. ACM SIGMOD Int. Conf. on Management  */
/* of Data, Philadelphia, PA, 1999."                                                */
/* -------------------------------------------------------------------------------- */

/*#include <stream.h>
#include <fstream.h>
#include <iomanip.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <errno.h>
#include <math.h>
#include <time.h>*/

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <cerrno>
#include <cmath>
#include <ctime>


using namespace std;

/* -------------------------------------------------------------------------------- */
/* Defines */
#define assignMin(a,b) if(a>b) a=b;
#define assignMax(a,b) if(a<b) a=b;
#define UNDEF -1
#define TRUE (1==1)
#define FALSE (1==0)
#define max(a,b) ((a)>(b)?(a):(b))
#define min(a,b) ((a)<(b)?(a):(b))
#define INFINITY (2*GD.eps)

// ======================================== 
// global datavalues
class globalData {
public:
  globalData(void) { dim=-1; }
  int dim;         // dimension of the data
  float eps;       // eps for the clustering algorithm
  int minPts;      // minPts value for the clustering algorithm
};

globalData GD;

// ======================================== 
// a data point
class point {
public:
  point(void) 
            { init(); }
  point(const point &p)
            { init(); assign(p); }
  ~point() 
            { delete[] f; }
  float getCoord(int _dim) const
            { assert((_dim>=0)&&(_dim<GD.dim)); return(f[_dim]); }
  void  setCoord(int _dim, float _f) 
            { assert((_dim>=0)&&(_dim<GD.dim)); f[_dim] = _f; }
  void  assign(const point &p) 
            { for(int i=0; i<GD.dim; ++i) f[i]=p.getCoord(i); }
  void  assign(float *_f) 
            { for(int i=0; i<GD.dim; ++i) f[i]=_f[i]; }
  void print(void) const
            { cout << "[ "; 
	      for(int i=0; i<GD.dim-1; ++i) {
		cout << f[i] << ", ";
	      }
	      cout << f[GD.dim-1] << "]\n" << flush;
	    }
  void fileread(FILE *infile)
	{ fread(f, sizeof(float), GD.dim, infile); assert(!feof(infile)); }
	
  float dist(const point &p) const
	
	
            { /* Euclidian Distance */
				float dist = 0, k = 0,kk = 0, minimum = 0,offset;
				
		for(int i=0; i<GD.dim; ++i) {
			if (f[i] < minimum) { minimum = f[i]; }
			if (p.f[i] < minimum) {minimum = p.f[i]; }
			if ((f[i] < -15) || (p.f[i] < -15)) { k = k + 1;}
			if ((f[i] < -15) && (p.f[i] < -15)) { kk = kk + 1;}
		}
										
			//better11 was -35, k = 0.1
			//better12 was -25, k = 0.2
		    //better13 was -25, k = 0.2 extra && comparison
		    //better15 was -25, k = 0.2 extra && comparison erase one offest from else
				
	      for(int i=0; i<GD.dim; ++i) {
			  //offset = abs( f[i] * ( f[i] - p.f[i]) ) ;
			
			 // if((f[i] == 0 && p.f[i] != 0)) {
			//	offset = (p.f[i]-f[i]);
			//	  offset = offset*offset;
			//  }
			//  else {
			//    offset = p.f[i]-f[i];
			//  }
			//  
			//  offset = offset*offset;
			  
			  
			  if ((f[i] < -15) && (p.f[i] < -15)) {
				  offset = (p.f[i]-f[i]);
				  offset = offset*offset;
				  offset = offset*offset * abs(minimum) / (k * (abs(minimum) + abs(min(f[i],p.f[i])))   ); 
				  
			  }
			  else if ((f[i] < -15) || (p.f[i] < -15)) { 
				  offset = (p.f[i]-f[i]);
				  offset = offset*offset;
				  offset = offset*offset* (k-kk)*(abs(minimum) +  abs(min(f[i],p.f[i]))) / (abs(minimum) * k*kk);
				  
		     //this is the original code... seems to be the best still
			 //		if ((f[i] < -25) && (p.f[i] < -25)) {
			 //			offset = (p.f[i]-f[i]);
			 // 			offset = offset*offset;
			 // 			offset = offset*offset/(kk* abs(min(f[i],p.f[i]))); 
			  
			 // 		}
			 // 		else if ((f[i] < -25) || (p.f[i] < -25)) { 
			//			  offset = (p.f[i]-f[i]);
			// 			  offset = offset*offset;
			//			  offset = offset*offset/(k * abs(min(f[i],p.f[i])));
		
			  
//		if ((f[i] < -25) && (p.f[i] < -25)) {
//			offset = (p.f[i]-f[i])/sqrt(kk* abs(min(f[i],p.f[i])));
//			offset = offset*offset;
//			offset = offset*offset; 
			
//		}
//		else if ((f[i] < -25) || (p.f[i] < -25)) { 
//			  offset = (p.f[i]-f[i])/sqrt(k * abs(min(f[i],p.f[i])));
//			  offset = offset*offset;
//			  offset = offset*offset;
			  //offset = offset*offset;
			  //offset = offset*offset;
		
		}
		else { offset = 0; 
		}
			  
			  
			  
			  
			  //if ((f[i] > -35) || (p.f[i] > -35)) {
				//  off = max(abs(f[i]),abs(p.f[i]));
				//  offset = offset * off * off * off* off / 35*35*35*35;
			  //}
			 // else { offset = 0; }
				  
			  
			  //offset = offset*offset;
			  //offset = offset*offset;
			  //offset = offset*offset;
			  
			  dist += offset;
	      }
	      return(sqrt(dist));
				
	    }
  float distSC(const point &p) const
            { // distSC(p1,p2) if distance Between p1 und p2 ist smaller 
	      // than or equal to GD.eps and INFINITY otherwise
	      float offset, eps;
				
				
	      //eps = GD.eps*GD.eps;
				
		//for(int i=0; i<GD.dim; ++i) {
		//offset = p.f[i]-f[i];
		//eps -= (offset*offset);
		//if(eps<0) return(INFINITY);
		//}
		//return(sqrt(GD.eps*GD.eps-eps));	
				
				
				//jan start
				eps = GD.eps;
				offset = dist(p); 
				if ((eps-offset)<0) return(INFINITY);
				return(offset);		
				//jan end
				
				
	     
	    }

private:
  inline void init(void) { f=new float[GD.dim]; }
  float  *f;     // values
  friend class mbr;
};

// ======================================== 
// a mbr (minimum bounding rectangle)     //I guess this is used only for lean scan
class mbr {
public:
  mbr(void) 
            { }
  mbr(const point &p) 
            { set(p); }
  ~mbr()
            { }
  void set(const point &p) 
            { setLB(p); setUB(p); }
  void setLB(const point &p) 
            { lb.assign(p); }
  void setUB(const point &p) 
            { ub.assign(p); }
  const point &getLB(void)
            { return(lb); }
  const point &getUB(void)
            { return(ub); }
  void include(const point &p) // resize the mbr so that it includes p
            { for(int i=0; i<GD.dim; ++i) {
                assignMin(lb.f[i], p.f[i]); 
		assignMax(ub.f[i], p.f[i]); 
	      } 
	    }
  float dist(const point &p)
            { float dist = 0, offset;
	      for(int i=0; i<GD.dim; ++i) {
		if(ub.f[i] <= p.f[i])
		  offset = p.f[i] - ub.f[i];
		else if(lb.f[i] >= p.f[i])
		  offset = lb.f[i] - p.f[i];
		else
		  offset = 0;
		dist += (offset*offset);
	      }
	      return(sqrt(dist));
	    }

private:
  point lb; // lower bound
  point ub; // upper bound
};

// ======================================== 
// a point list
#define PTLISTLEN 20
class ptListContent {
public:
  int pts; // distance of point pts 
  float d; // is d
};

class ptListEntry {
public:
  ptListEntry(void) { next = NULL; }
  ~ptListEntry()    { if(next) delete next; }
private:
  friend class ptList;
  friend class ptListIter;
  ptListContent pc[PTLISTLEN];
  ptListEntry *next;
};

class ptList {
public:
  ptList(void)   { nextfree = 0; entry = NULL; }
  ~ptList()      { if(entry) delete entry; }
  void add(ptListContent *pc)
            { if(entry==NULL) { // list is completely empty
		entry = new ptListEntry;
	      } else if(nextfree==PTLISTLEN) { // current entry is full
		ptListEntry *newentry = new ptListEntry;
		newentry->next = entry;
		entry = newentry;
		nextfree = 0;
	      }
	      entry->pc[nextfree] = *pc;
	      ++nextfree;
	    }

private:
  friend class ptListIter;
  int nextfree;
  ptListEntry *entry;
};

class ptListIter {
public:
  ptListIter(const ptList *pl)
            { maxEntry = pl->nextfree; entry = pl->entry; i=0; }
  ptListContent *next(void) 
            { if(entry==NULL) return(NULL);
  	      if(i<maxEntry) {
		++i;
		return(&(entry->pc[i-1]));
	      } else {
		entry = entry->next;
		i=0;
		maxEntry = PTLISTLEN;
		return(next());
	      }
	    }
private:
  int maxEntry, i;
  ptListEntry *entry;
};

// ======================================== 
// an entry in a datapage
class dataBaseEntry {
public:
  dataBaseEntry()              { }
  ~dataBaseEntry()             { }
  void setPoint(point &_p)     { p.assign(_p);  }
  const point &getPoint(void)  { return(p); }
private:
  point p;        // the data point
};

// ======================================== 
// a datapage
class dataPage {
public:
  dataPage(void)
            { numEntries = 0; entry = NULL; }
  ~dataPage()
            { if(entry!=NULL) delete[] entry; }
  void setNumEntries(int _numEntries)
            { numEntries = _numEntries;
	      if(entry!=NULL) delete[] entry;
	      entry = new dataBaseEntry[numEntries]; 
	    }
  int             firstPointIndex;   // number of the first point in this page
  int             numEntries;        // how many entries does this data page have
  dataBaseEntry*  entry;             // array of size "numEntries"
  mbr             rect;              // minimal bounding rectangle of this page
private:
};

// ======================================== 
// the database - abstract base class
class dataBase {
public:
  dataBase(void) { }
  virtual ~dataBase() { }
  int                      getNumPoints(void) { return numPoints; }
  virtual const point &    getPoint(int index) = 0;
  virtual ptList *         rangeQuery(const point &p) = 0;
  virtual int              specialNNQuery(int *handled) = 0;
protected:
  int numPoints;
};

// ======================================== 
// a very simple database, no index, just the points
class simpleDataBase : public dataBase {
public:
  simpleDataBase(void);
  ~simpleDataBase();
  int              load(const char *filename, int _numPoints);
  const point &    getPoint(int index) { assert(index<numPoints); return(entry[index].getPoint()); }
  ptList *         rangeQuery(const point &p);
  int              specialNNQuery(int *handled);
private:
  dataBaseEntry *entry;
};

// ======================================== 
// a database using a lean-tree
class leanTreeDataBase : public dataBase {
public:
  leanTreeDataBase(void);
  ~leanTreeDataBase();
  int load(const char *filename);
  const point &getPoint(int index)
            { int p = pageContainingPoint(index);
	      return(page[p].entry[index-page[p].firstPointIndex].getPoint());
	    }
  ptList *         rangeQuery(const point &p);
  int              specialNNQuery(int *handled);
  int              getNumPages(void) { return numPages; }
  int              getPageSize(void) { return pageSize; }
  int              getMaxEntriesPerPage(void) { return maxEntriesPerPage; }
private:
  int pageContainingPoint(int index) // compute pageNum containing point index
            { assert(index<numPoints);
	      int p=0, s=0;
	      while(s+page[p].numEntries-1<index) { s += page[p].numEntries; ++p; }
	      return(p);
	    }
  int numPages;
  int maxEntriesPerPage;
  int pageSize;
  dataPage *page;
};

// ======================================== 
// the OPTICS algorithm
#define HANDLED -1
#define NOTHANDLED -2
// ======================================== 
// helper class to manage the seed list
class seedList {
public:
  seedList(int _size);
  ~seedList();
  void addRandomPoint(dataBase *db);
  void update       (int index, float r);
  void getMin       (int *index, float *r);
  int  getNumEntries(void)  { return(numEntries); }
private:
  void adjustHeap   (int index);
  int     numEntries;
  int     doneIndex;
  int     size;
  float * slReach;
  int *   slIndex;
  int *   slInvIndex;
};

// the optics algorithm
class optics {
public:
  optics               (dataBase *_db, ostream *_os);
  ~optics              ();
  float updateAllReach (int index);
  void output          (int index, float coreLevel, float reach);
private:
  dataBase *db;
  ostream *os;
  seedList *sl;
  float *ndist;
  char buffer[256];
};

// ================================================================================
// ========== main                                                       ==========
// ================================================================================
int main(int argc, char **argv) {
  int runType;
  clock_t Tstart, Tend;

  printf("\n");
  printf("=====================================================================\n");
  printf("| OPTICS                                                            |\n");
  printf("| (c) University of Munich, Database Group                          |\n");
  printf("=====================================================================\n\n");

  if(argc<2) {
    printf("**** Usage: %s -s <eps> <minPts> <dim> <n> <file> <outfile>\n", argv[0]);
    printf("            Use sequential scan (.dump file) with\n");
    printf("            <eps>      = epsilon value for clustering algorithm\n"); 
    printf("            <minPts>   = minPts value for clustering algorithm\n");
    printf("            <pageSize> = page size\n");
    printf("            <dim>      = dimension of the data set\n");
    printf("            <n>        = number of points in the data set\n");
    printf("            <file>     = filename of the data set (.dump)\n");
    printf("            <outfile>  = filename of the output file (.cri)\n");
    printf("     <<OR>>\n");
    printf("**** Usage: %s -l <eps> <minPts> <file> <outfile>\n", argv[0]);
    printf("            Use lean-tree (.ltree file) with\n");
    printf("            <eps>      = epsilon value for clustering algorithm\n"); 
    printf("            <minPts>   = minPts value for clustering algorithm\n");
    printf("            <file>     = filename of the lean-tree file (.ltree)\n");
    printf("            <outfile>  = filename of the output file (.cri)\n");
    printf("     <<OR>>\n");
    printf("\n");
    exit(-1);
  }
  if(argv[1][0]=='-') {
    runType = argv[1][1];
  } else {
    printf("=== FATAL ERROR ===: wrong param 1. Run without params for usage info.\n");
    exit(-3);
  }

  if(runType=='s') {
    /* use dump-file as input */
    if(argc!=8) {
      cout << "=== FATAL ERROR ===: wrong params. Run without params for usage info.\n";
      exit(-3);
    }

    GD.eps         = atof(argv[2]);
    GD.minPts      = atoi(argv[3]);
    GD.dim         = atoi(argv[4]);
    int numPoints  = atoi(argv[5]);
    char *datafile = argv[6];
    char *outfile  = argv[7];

    simpleDataBase db;
    if(db.load(datafile, numPoints)!=0) exit(-2);
    printf("\n**** Using the following parameters:\n");
    printf("     eps         = %f\n",GD.eps);
    printf("     minPts      = %d\n",GD.minPts);
    printf("     dim         = %d\n",GD.dim);
    printf("     n           = %d\n",db.getNumPoints());
    ofstream os(outfile);                                          assert(os);
    
    printf("\n**** Executing OPTICS (based on sequential scan)\n");
    Tstart = clock();
    optics *op = new optics(&db, &os);                   assert(op);
    Tend = clock();
    delete op;
  
    int CPUtime= (int) ((Tend-Tstart)/CLOCKS_PER_SEC);
    printf("\n\n**** Performance\n");
    printf("     total cpu time used      : %d sec.\n", CPUtime);
    
    os.close();
  } else if(runType=='l') {
    /* use lean tree as input */
    if(argc!=6) {
      cout << "=== FATAL ERROR ===: wrong params. Run without params for usage info.\n";
      exit(-3);
    }

    GD.eps         = atof(argv[2]);
    GD.minPts      = atoi(argv[3]);
    char *datafile = argv[4];
    char *outfile  = argv[5];

    leanTreeDataBase db;
    if(db.load(datafile)!=0) exit(-2);
    printf("\n**** Using the following parameters:\n");
    printf("     eps         = %f\n",GD.eps);
    printf("     minPts      = %d\n",GD.minPts);
    printf("     dim         = %d\n",GD.dim);
    printf("     n           = %d\n",db.getNumPoints());
    printf("     PageSize    = %d\n",db.getPageSize());
    printf("       EntriesPerPage (max) = %d\n", db.getMaxEntriesPerPage());
    printf("       numDataPages         = %d\n", db.getNumPages());
    ofstream os(outfile); assert(os);
    
    printf("\n**** Executing OPTICS (based on lean tree)\n");
    Tstart = clock();
    optics *op = new optics(&db, &os); assert(op);
    Tend = clock();
    delete op;
    
    int CPUtime= (int) ((Tend-Tstart)/CLOCKS_PER_SEC);
    printf("\n\n**** Performance\n");
    printf("     total cpu time used      : %d sec.\n", CPUtime);
    
    os.close();
  } else {
    cout << "=== FATAL ERROR ===: unknown param 1. Run without params for usage info.\n";
    exit(-3);
  }

  return(0);
}


// ================================================================================
// ========== class optics                                               ==========
// ================================================================================

// --------------------------------------------------------------------------------
// execute the optics algorithm
optics::optics(dataBase *_db, ostream *_os) {
  int donePoints = 0;
  int index;
  float r, coreLevel;

  db = _db;
  os = _os;
  sl = new seedList(db->getNumPoints()); assert(sl);
  ndist = new float[GD.minPts]; assert(ndist);
  while(donePoints<db->getNumPoints()) {
    if(sl->getNumEntries()==0) {
      // pick a random point and add to seedlist with r=infinity
      sl->addRandomPoint(db);
    }
    sl->getMin(&index, &r); // get and delete point with smallest r from seedlist
    coreLevel = updateAllReach(index); // update seedlist with reachabilities from point index
    output(index, coreLevel, r); // output this point with it's current reachability
    ++donePoints;
	  if((donePoints%100) ==0) {																//Jan
		  cout << "\r\t"<<(int)(donePoints*100.0/db->getNumPoints()) << "% [" 
		  << donePoints << " of " << db->getNumPoints() << "]" << flush;
	  }
  }
}

// --------------------------------------------------------------------------------
// free memory
optics::~optics() {
  delete sl;
  delete[] ndist;
}

// --------------------------------------------------------------------------------
// enter/update all the reachabilities of all points as seen from point index
// standard version of OPTICS (non-join based)
float optics::updateAllReach(int index) {
  ptListContent *pc;
  int i;

  for(i=0; i<GD.minPts; ++i) 
    ndist[i] = INFINITY;

  // execute range query
  ptList *pl = db->rangeQuery(db->getPoint(index));

  // compute core level
  ptListIter ptli(pl);
  while((pc = ptli.next())) {
    if(ndist[GD.minPts-1]>pc->d) {
      // have to include this dist into ndist array
      i = GD.minPts-2;
      while((ndist[i]>pc->d)&&(i>=0)) {
	ndist[i+1] = ndist[i];
	--i;
      }
      ndist[i+1] = pc->d;
    }
  }
  float coreLevel = ndist[GD.minPts-1];

  // update all points in eps-neighborhood
  if(coreLevel<INFINITY) {
    // point is core point
    ptListIter ptli(pl);
    while((pc = ptli.next())) {
      sl->update(pc->pts, max(pc->d, coreLevel));
    }
  }
  delete pl;
  return(coreLevel);
}

// --------------------------------------------------------------------------------
// output the point index with reachability r 
void optics::output(int index, float c, float r) {
  // cannot use stream functions to print floats, because
  // even if we set ios::fixed, it sometimes prints e-notation.
  const point &p = db->getPoint(index);

  for(int i=0; i<GD.dim; ++i) {
    sprintf(buffer, "%10.8f ", p.getCoord(i));
    (*os) << buffer;
  }
  sprintf(buffer, "%10.8f %10.8f \"%d\"", c, ((r==INFINITY)?-1:r), index);
  (*os) << buffer << endl;
}

// ================================================================================
// ========== class seedList                                             ==========
// ================================================================================

// init all data structures
seedList::seedList(int _size) {
  size       = _size;
  numEntries = 0;
  doneIndex  = 0;
  slReach    = new float[size];   assert(slReach);
  slIndex    = new int[size];     assert(slIndex);
  slInvIndex = new int[size];     assert(slInvIndex);
  for(int i=0; i<size; ++i)
    slIndex[i] = NOTHANDLED;
}

seedList::~seedList() {
  delete[] slReach;
  delete[] slIndex;
  delete[] slInvIndex;
}

void seedList::addRandomPoint(dataBase *db) {
  // do NN-query around origin
  int index = db->specialNNQuery(slIndex);
  update(index, INFINITY);
}

// adjust the heap structure, assuming only position "index" IN THE HEAP has been
// updated/added
void seedList::adjustHeap(int index) {
  int j;
  float hf;
  int hi;

  // do we need to bubble up or sink down?
  if(index>0) { // we are not the root
    j = (index-1)/2; // j is father
    if(slReach[index]<slReach[j]) {
      // bubble up
      while(index>0) {
	// exchange(index, j);
	hf = slReach[index];
	slReach[index] = slReach[j];
	slReach[j] = hf;

	slIndex[slInvIndex[index]] = j;
	slIndex[slInvIndex[j]] = index;

	hi = slInvIndex[index];
	slInvIndex[index] = slInvIndex[j];
	slInvIndex[j] = hi;
	// end exchange(index, j);
	index = j;
	if(index>0) {
	  j = (index-1)/2;
	  if(slReach[index]>=slReach[j]) // heap ok
	    index = 0;
	}
      }
      return; // STOP HERE
    }
  }
  // sink down if necessary
  j = 2*index+1; // left son
  while(j < numEntries) { // we are not a leaf
    if(j<numEntries-1)  // right leaf does exist
      if(slReach[j+1]<slReach[j]) // right leaf is smaller
	++j;
    // slReach[j] is smaller son if slReach[index]!
    if(slReach[j]<slReach[index]) { // heap not ok
      // exchange(index, j);
      hf = slReach[index];
      slReach[index] = slReach[j];
      slReach[j] = hf;
      
      slIndex[slInvIndex[index]] = j;
      slIndex[slInvIndex[j]] = index;

      hi = slInvIndex[index];
      slInvIndex[index] = slInvIndex[j];
      slInvIndex[j] = hi;
      // end exchange(index, j);
      index = j;
      j = index*2+1;
    } else  // heap ok
      j = numEntries;
  }
}

// add or update point "index" with reachability value "r"
void seedList::update(int index, float r) {
  if(slIndex[index]!=HANDLED) {
    if(slIndex[index]==NOTHANDLED) {
      // point needs to be added to the seedlist
      slIndex[index] = numEntries;
      slReach[numEntries] = r;
      slInvIndex[numEntries] = index;
      numEntries++;
      adjustHeap(numEntries-1);
    } else {
      // point is already in the seedlist
      if(slReach[slIndex[index]] > r) {
	// need to update because new reachability is smaller
	slReach[slIndex[index]] = r;
	adjustHeap(slIndex[index]);
      }
    }
  }
}

// get and delete the point with the smallest r-value from the seedlist
void seedList::getMin(int *index, float *r) {
  *r = slReach[0];
  *index = slInvIndex[0];
  --numEntries;
  slIndex[slInvIndex[0]] = HANDLED;
  if(numEntries>0) {
    // do this only if the list is not empty
    slReach[0] = slReach[numEntries];
    slInvIndex[0] = slInvIndex[numEntries];
    slIndex[slInvIndex[0]] = 0;
    adjustHeap(0);
  }
}


// ================================================================================
// ========== class simpleDataBase                                             ==========
// ================================================================================

// --------------------------------------------------------------------------------
// prepare database to be read into memory
simpleDataBase::simpleDataBase(void) {
  entry = NULL;
}

// --------------------------------------------------------------------------------
// if a database was read into memory, free the memory
simpleDataBase::~simpleDataBase() {
  if(entry!=NULL) {
    delete[] entry;
  }
}

// --------------------------------------------------------------------------------
// load the binary data in "datafile". return 0 if ok, something else otherwise
int simpleDataBase::load(const char *filename, int _numPoints) {
  FILE *infile;
  int i;
  point *p;

  printf("**** Reading Input File (.dump)\n");
  numPoints = _numPoints;
  infile = fopen(filename, "rb");                           assert(infile);
  p = new point;                                            assert(p);
  entry = new dataBaseEntry[numPoints];                     assert(entry);

  for(i=0; i<numPoints; ++i) {  
    if((i%1000) ==0) cout << "\r " << i << flush;
    p->fileread(infile);
    entry[i].setPoint(*p);
  }
  printf("\r     #points read = %d\n",numPoints);
  delete p;
  return(0);
}

// --------------------------------------------------------------------------------
// execute a range query on the database
// NOTE: the returned ptList has to be deleted by the caller!
ptList * simpleDataBase::rangeQuery(const point &p) {
  ptList *pl = new ptList;
  ptListContent pc;

  for(int i=0; i<numPoints; ++i) { // for all points
    if((pc.d = p.distSC(entry[i].getPoint()))<INFINITY) {
      pc.pts = i;
      pl->add(&pc);
    }
  }
  return(pl);
}

// --------------------------------------------------------------------------------
// do a special NN-query: compute point next to origin that has not been
// handled
int simpleDataBase::specialNNQuery(int *handled) {
  float currdist=0, newdist;
  int currindex = -1, i;
  point origin;

  for(i=0; i<GD.dim; ++i) origin.setCoord(i, 0);

  for(i=0; i<numPoints; ++i) { // for evert point
    if(handled[i]==NOTHANDLED) {
      newdist = origin.dist(entry[i].getPoint());
      if((newdist<currdist)||(currindex==-1)) {
	currdist = newdist;
	currindex = i;
      }
    }
  }
  return(currindex);
}



// ================================================================================
// ========== class leanTreeDataBase                                             ==========
// ================================================================================

// --------------------------------------------------------------------------------
// prepare database to be read into memory
leanTreeDataBase::leanTreeDataBase(void) {
  page = NULL;
}

// --------------------------------------------------------------------------------
// if a database was read into memory, free the memory
leanTreeDataBase::~leanTreeDataBase() {
  if(page!=NULL) {
    delete[] page;
  }
}

// --------------------------------------------------------------------------------
// load the lean tree data in "datafile". return 0 if ok, something else otherwise
int leanTreeDataBase::load(const char *filename) {
  FILE *infile;
  int i,j, nRead=0, numEntries;
  point *p;

  printf("**** Reading Input File (.ltree)\n");
  infile = fopen(filename, "r");                             assert(infile);
  /* read dimension */
  fread(&GD.dim, sizeof(int), 1, infile);                    assert(!feof(infile));
  /* read entriesPerPage */
  fread(&maxEntriesPerPage, sizeof(int), 1, infile);         assert(!feof(infile));
  /* read number of data points */
  fread(&numPoints, sizeof(int), 1, infile);                 assert(!feof(infile));
  /* read numDataPages */
  fread(&numPages, sizeof(int), 1, infile);                  assert(!feof(infile));

  pageSize = (maxEntriesPerPage+2)*GD.dim*sizeof(float)+sizeof(int);

  p = new point;
  page = new dataPage[numPages];                             assert(page);
  // read the data pages
  for(i=0; i<numPages; ++i) {
    // read numEntries for this page
    fread(&numEntries, sizeof(int), 1, infile);              assert(!feof(infile));
    page[i].setNumEntries(numEntries);
    page[i].firstPointIndex = nRead;
    // read lb and ub (mbr) from file
    p->fileread(infile);  page[i].rect.setLB(*p);
    p->fileread(infile);  page[i].rect.setUB(*p);
    // read the points
    for(j=0; j<numEntries; ++j) {
      if((nRead%1000)==0) cout << "\rpages=" << i << " points=" << nRead << flush;
      ++nRead;
      p->fileread(infile);
      page[i].entry[j].setPoint(*p);
      page[i].rect.include(*p);
    }
  }
  printf("\r     numPagesRead      = %d\n",numPages);
  printf("     numDatapointsRead = %d\n",nRead);
  delete p;
  return(0);
}

// --------------------------------------------------------------------------------
// execute a range query on the database
// NOTE: the returned ptList has to be deleted by the caller!
ptList * leanTreeDataBase::rangeQuery(const point &p) {
  ptList *pl = new ptList;
  ptListContent pc;

  for(int i=0; i<numPages; ++i) { // for all pages
    if(page[i].rect.dist(p)<=GD.eps) { // page may contain points in range
      for(int j=0; j<page[i].numEntries; ++j) {
	if((pc.d=p.distSC(page[i].entry[j].getPoint()))<INFINITY) {
	  pc.pts = page[i].firstPointIndex+j;
	  pl->add(&pc);
	}
      }
    }
  }
  return(pl);
}

// --------------------------------------------------------------------------------
// do a special NN-query: compute point next to origin that has not been
// handled
int leanTreeDataBase::specialNNQuery(int *handled) {
  float currdist=0, newdist;
  int index = -1, i, j;
  point origin;

  for(i=0; i<GD.dim; ++i) origin.setCoord(i, 0);

  for(i=0; i<numPages; ++i) { 
    // for every page
    if((page[i].rect.dist(origin) < currdist)||(index==-1)) {
      // consider page
      for(j=0; j<page[i].numEntries; ++j) {
	// for every point in the page
	if(handled[page[i].firstPointIndex+j]==NOTHANDLED) {
	  // point has not already been handled
	  newdist = origin.dist(page[i].entry[j].getPoint());
	  if((newdist<currdist)||(index==-1)) {
	    currdist = newdist;
	    index = page[i].firstPointIndex + j;
	  }
	}
      }
    }
  }
  return(index);
}


// ================================================================================
// ========== end of file                                                ==========
// ================================================================================
