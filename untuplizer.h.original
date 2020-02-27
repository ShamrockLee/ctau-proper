/* Reader of ggNtuplizer's and other TTrees.
 *
 * Works both with old as well as with new ntuples. In particular, for tree
 * branches containing vector<...> objects of elementary data types as well as
 * two-dimensional vector<vector<float> > or vector<vector<int> > arrays,
 * underlying arrays are returned (vector::front()) by TreeReader::Get*() call
 * family.
 *
 * NOTE: vector<bool> cannot be handled in this way since booleans are packed
 * as bits into bytes.
 *
 * NOTE: for branches containing variable length arrays, corresponding counter
 * branches must be read off first.
 *
 * Usage of int, long, ... elementary data types for retrieving branch contents
 * is not portable and should be avoided. This is due to the fact that
 * sizeof(int) and sizeof(long) may differ on 32bit and 64bit architectures.
 * Therefore, it is advised to use the following conventions for types of branch
 * contents:
 *   1. Int_t instead of int;
 *   2. UInt_t instead of unsigned int;
 *   3. Long64_t instead of long;
 *   4. ULong64_t instead of unsigned long, etc.
 *
 * For instance, it is suggested to use
 *   Int_t nEle = reader.GetInt("nEle");
 * instead of
 *   int nEle = reader.GetInt("nEle");
 *
 * It is generally safe, however, to use float and double in place of Float_t
 * and Double_t.
 *
 * Some notes:
 * 1. In-place modifications of arrays (within array boundaries) are allowed and
 *    encouraged.
 * 2. The code loads only requested tree branches. Therefore, in general, an
 *    analysis will run much faster.
 * 3. Branch contents are accessed directly by accessing leafs. Thus, it is
 *    assumed that all leafs have unique names (equal to corresponding names of
 *    branches) which is indeed the case when each branch contains only one
 *    leaf.
 * 4. In the code, no protection against invalid code usage is provided.
 * 5. Tree friends are not supported.
 */

#ifndef UNTUPLIZER_H
#define UNTUPLIZER_H

#include <map>
#include <string>
#include <vector>
#include <TFile.h>
#include <TChain.h>
#include <TLeafF.h>
#include <TLeafD.h>
#include <TLeafB.h>
#include <TLeafS.h>
#include <TLeafI.h>
#include <TLeafL.h>
#include <TLeafO.h>
#include <TSystem.h>
#include <TLeafObject.h>
#include <TLeafElement.h>

// prints a message and exits gracefully
#ifndef FATAL
#define FATAL(msg) do { fprintf(stderr, "FATAL: %s\n", msg); gSystem->Exit(1); } while (0)
#endif

class TreeReader {

 public:

   // types of branch/leaf contents
   // NOTE: some unsigned types may not be supported
   enum ETypes {
      kBool,            // single 1-byte boolean (TLeafO)
      kChar,            // single 1-byte integer (TLeafB)
      kShort,           // single 2-byte integer (TLeafS)
      kInt,             // single 4-byte integer (TLeafI)
      kFloat,           // single 4-byte float (TLeafF)
      kDouble,          // single 8-byte float (TLeafD)
      kLong64,          // single 8-byte integer (TLeafL)
      kArrBool,         // array of 1-byte booleans (TLeafO only)
      kArrChar,         // array of 1-byte integers (TLeafB or vector<char>)
      kArrCharTLeaf,    // array of 1-byte integers (TLeafB)
      kArrCharVector,   // array of 1-byte integers (vector<char>); NOTE: char=signed char is assumed
      kArrUCharVector,  // array of 1-byte unsigned integers (vector<unsigned char>)
      kArrShort,        // array of 2-byte integers (TLeafS or vector<short>)
      kArrShortTLeaf,   // array of 2-byte integers (TLeafS)
      kArrShortVector,  // array of 2-byte integers (vector<short>)
      kArrUShortVector, // array of 2-byte unsigned integers (vector<unsigned short>)
      kArrInt,          // array of 4-byte integers (TLeafI or vector<int>)
      kArrIntTLeaf,     // array of 4-byte integers (TLeafI)
      kArrIntVector,    // array of 4-byte integers (vector<int>)
      kArrUIntVector,   // array of 4-byte unsigned integers (vector<unsigned int>)
      kArrFloat,        // array of 4-byte floats (TLeafF or vector<float>)
      kArrFloatTLeaf,   // array of 4-byte floats (TLeafF)
      kArrFloatVector,  // array of 4-byte floats (vector<float>)
      kArrLong64,       // array of 8-byte integers (TLeafL or vector<long>)
      kArrLong64TLeaf,  // array of 8-byte integers (TLeafL)
      kArrLong64Vector, // array of 8-byte integers (vector<long>)
      kArrULong64Vector,// array of 8-byte unsigned integers (vector<unsigned long>)
      kArrVectInt,      // array of vector<int> (vector<vector<int> >)
      kArrVectFloat,    // array of vector<float> (vector<vector<float> >)
      kArrString,       // array of string 
      kArrStringVector, // array of string (vector<string>)
      kTObject,         // general object inherited from TObject
      kVoidPtr          // all other data types
   };

   // TTree/TChain initializers
   TreeReader(TTree* tree);
   TreeReader(const char* path, const char* treename = "tree/treeMaker");
   TreeReader(const char** paths, int npaths, const char* treename = "tree/treeMaker");
   TreeReader(std::vector<std::string> paths, const char* treename = "tree/treeMaker");

   virtual ~TreeReader();

   // simple getters
   TTree*   GetTree()        { return fTree; }
   Long64_t GetEntriesFast() { return fTree->GetEntriesFast(); }
   Bool_t   HasMC()          { return fkMC; }

   // useful to determine which type of variable to use for which branch
   void Print();

   // sets event number to retrieve next time TreeReader::Get*() called
   void GetEntry(Long64_t entry);

   // returns either a pointer to first element of one- or multi-dimensional
   // array or a pointer to an unprocessed object.
   void* GetPtr(const char* branch_name, ETypes cktype = kVoidPtr, Int_t* nsize = NULL);

   // specializations of the call above;
   // NOTE: use same methods with type casting for unsigned data types
   Char_t*    GetPtrChar   (const char* bname) { return (Char_t*)    GetPtr(bname, kArrChar);    }
   Short_t*   GetPtrShort  (const char* bname) { return (Short_t*)   GetPtr(bname, kArrShort);   }
   Int_t*     GetPtrInt    (const char* bname) { return (Int_t*)     GetPtr(bname, kArrInt);     }
   Float_t*   GetPtrFloat  (const char* bname) { return (Float_t*)   GetPtr(bname, kArrFloat);   }
   Long64_t*  GetPtrLong64 (const char* bname) { return (Long64_t*)  GetPtr(bname, kArrLong64);  }
   TObject*   GetPtrTObject(const char* bname) { return (TObject*)   GetPtr(bname, kTObject);    }

   // NOTE: this works only for old ntuples
   Bool_t*    GetPtrBool(const char* bname) { return (Bool_t*) GetPtr(bname, kArrBool); }

   // return branch values for elementary types
   Bool_t   GetBool  (const char* bname) { return ((Bool_t*)   GetPtr(bname, kBool))  [0]; }
   Char_t   GetChar  (const char* bname) { return ((Char_t*)   GetPtr(bname, kChar))  [0]; }
   Short_t  GetShort (const char* bname) { return ((Short_t*)  GetPtr(bname, kShort)) [0]; }
   Int_t    GetInt   (const char* bname) { return ((Int_t*)    GetPtr(bname, kInt))   [0]; }
   Float_t  GetFloat (const char* bname) { return ((Float_t*)  GetPtr(bname, kFloat)) [0]; }
   Double_t GetDouble(const char* bname) { return ((Double_t*) GetPtr(bname, kDouble))[0]; }
   Long64_t GetLong64(const char* bname) { return ((Long64_t*) GetPtr(bname, kLong64))[0]; }

   std::string* GetPtrString(const char* bname) { return (std::string*) GetPtr(bname, kArrString); }
   Int_t GetPtrStringSize(){return fStringVectorSize;}


   // vector<vector<float> > and vector<vector<int> > tree branches
   std::vector<Float_t>* GetPtrVectorFloat(const char* bname, Int_t &nsize) {
      return (std::vector<Float_t>*) GetPtr(bname, kArrVectFloat, &nsize);
   }
   std::vector<Float_t>* GetPtrVectorFloat(const char* bname) {
      return (std::vector<Float_t>*) GetPtr(bname, kArrVectFloat);
   }
   std::vector<Int_t>* GetPtrVectorInt(const char* bname, Int_t &nsize) {
      return (std::vector<Int_t>*) GetPtr(bname, kArrVectInt, &nsize);
   }
   std::vector<Int_t>* GetPtrVectorInt(const char* bname) {
      return (std::vector<Int_t>*) GetPtr(bname, kArrVectInt);
   }

 protected:

   void  InitSingleTTree(const char* path, const char* treename);
   void  InitTChain(const char** paths, int npaths, const char* treename);
   void  FindLeaf(const char* bname);

   TFile*   fFile;     // file handle associated with fTree
   TTree*   fTree;     // reference to TTree or TChain to read entries from
   Int_t    fTreeNum;  // (for TChain) current tree number in list of TTrees
   Bool_t   fkMC;      // if kTRUE, MC truth info is available
   Long64_t fEntry;    // fTree entry number to read with Get*() methods
   Int_t    fStringVectorSize; 

   // caching
   std::map<std::string,int> fLeafIdx; // leaf name => index in arrays below
   std::vector<TLeaf*>  fLeafAddr;     // cached leaf address
   std::vector<ETypes>  fLeafType;     // cached type of leaf contents
   std::vector<void*>   fLeafValue;    // cached payload address
};

//______________________________________________________________________________
TreeReader::TreeReader(TTree* tree) :
   fFile(0),
   fTree(0),
   fTreeNum(-1),
   fkMC(kFALSE),
   fEntry(-1),
   fStringVectorSize(0)
{
   /* Associates an external TTree or TChain with this class.
    *
    * tree = reference to TTree or TChain. The caller is the owner of the
    * object.
    */

   // verify sizes of elementary data types
   if (sizeof(int) != 4 || sizeof(long) != 8 || sizeof(float) != 4)
      FATAL("int/long/float of unsupported size");

   fTree = tree;

   // find out availability of MC truth info (check existence of "nMC" branch)
   fkMC = fTree->GetBranch("nMC") ? kTRUE : kFALSE;
}

//______________________________________________________________________________
TreeReader::TreeReader(const char* path, const char* treename) :
   fFile(0),
   fTree(0),
   fTreeNum(-1),
   fkMC(kFALSE),
   fEntry(-1)
{
   /* Takes TTree from a root file.
    *
    * path = any root-supported path to a file with TTree.
    */

   InitSingleTTree(path, treename);
}

//______________________________________________________________________________
TreeReader::TreeReader(const char** paths, int npaths, const char* treename) :
   fFile(0),
   fTree(0),
   fTreeNum(-1),
   fkMC(kFALSE),
   fEntry(-1)
{
   /* Makes TChain of trees from several root files.
    *
    * paths = array of any root-supported paths to root files with TTrees;
    * npaths = number of elements in the array above.
    */

   InitTChain(paths, npaths, treename);
}

//______________________________________________________________________________
TreeReader::TreeReader(std::vector<std::string> paths, const char* treename) :
   fFile(0),
   fTree(0),
   fTreeNum(-1),
   fkMC(kFALSE),
   fEntry(-1)
{
   /* Makes TChain of trees from several root files.
    *
    * paths = array of any root-supported paths to root files with TTrees.
    */

   int npaths = (int) paths.size();
   const char* paths2[npaths];

   // convert vector<string> into array of null-terminated strings
   for (int i = 0; i < npaths; i++)
       paths2[i] = paths[i].c_str();

   InitTChain(paths2, npaths, treename);
}

//______________________________________________________________________________
TreeReader::~TreeReader()
{
   /* Frees memory.
    */

   if (fTree) delete fTree;
   if (fFile) delete fFile;
}

//______________________________________________________________________________
void TreeReader::Print()
{
   /* Prints to stdout branch names together with descriptions of their content
    * as to be used in an analysis code.
    *
    * Useful to determine which type of variable to use for which branch.
    *
    * NOTE: for branches containing vector<...> objects of elementary data types
    * or vector<...> objects of vector<int> or vector<float> data types, types
    * of underlying arrays are shown (as to be returned by TreeReader::Get*()
    * call family).
    */

   Printf("Branch names with content descriptions:");

   // access leafs directly
   TObjArray* leafs = fTree->GetListOfLeaves();

   // leaf loop
   for (int i = 0; i < leafs->GetEntriesFast(); i++) {
      TLeaf* leaf = dynamic_cast<TLeaf*>(leafs->At(i));
      if (!leaf)
         FATAL("leaf is NULL");

      // vector<...> and other objects not inherited from TObject
      if (leaf->IsA() == TLeafElement::Class()) {
         // content descriptor
         std::string descr(leaf->GetBranch()->GetClassName());

         // known one- and two-dimensional vector<...> arrays
         // NOTE: force showing data types of well-defined sizes
         const char* types[] = {
            "float", "Float_t", "char", "Char_t", "signed char", "unsigned char", "UChar_t",
            "short", "short int", "signed short", "signed short int", "Short_t",
            "unsigned short", "unsigned short int", "UShort_t",
            "int", "signed", "signed int", "Int_t", "unsigned int", "unsigned", "UInt_t",
            "long", "long int", "signed long", "signed long int", "Long64_t", "Long_t",
            "unsigned long", "unsigned long int", "ULong64_t", "ULong_t",
            "vector<float> ", "vector<int> ","vector<std::string> "};
         const char* type_descr[] = {
            "float", "float", "Char_t", "Char_t", "Char_t", "UChar_t", "UChar_t",
            "Short_t", "Short_t", "Short_t", "Short_t", "Short_t",
            "UShort_t", "UShort_t", "UShort_t",
            "Int_t", "Int_t", "Int_t", "Int_t", "UInt_t", "UInt_t", "UInt_t",
            "Long64_t", "Long64_t", "Long64_t", "Long64_t", "Long64_t", "Long64_t",
            "ULong64_t", "ULong64_t", "ULong64_t", "Long64_t",
            "vector<float>", "vector<Int_t>", "vector<std::string>"};

         for (int c = 0; c < 34; c++)
            if (descr.compare(Form("vector<%s>", types[c])) == 0) {
               descr = Form("%s*", type_descr[c]);
               break;
            }

         Printf("  %-36s: %s", leaf->GetName(), descr.data());
      }

      // objects inherited from TObject
      else if (leaf->IsA() == TLeafObject::Class()) {
         // content descriptor
         std::string descr(leaf->GetBranch()->GetClassName());

         Printf("  %-36s: %s", leaf->GetName(), descr.data());
      }

      // elementary data types; fixed/variable length arrays of elementary types
      else {
         std::string descr;

         // NOTE: only necessary leaf types are listed
         if (leaf->IsA() == TLeafF::Class())
            descr = "Float_t";
         else if (leaf->IsA() == TLeafI::Class())
            descr = "Int_t";
         else if (leaf->IsA() == TLeafL::Class())
            descr = "Long64_t";
         else if (leaf->IsA() == TLeafO::Class())
            descr = "Bool_t";
         else if (leaf->IsA() == TLeafD::Class())
            descr = "Double_t";
         else if (leaf->IsA() == TLeafB::Class())
            descr = "Char_t";
         else if (leaf->IsA() == TLeafS::Class())
            descr = "Short_t";
         else
            FATAL(Form("unsupported leaf of class %s", leaf->ClassName()));

         // single elementary variable vs fixed/variable length array
         if (!leaf->GetLeafCount() && leaf->GetLenStatic() == 1)
            Printf("  %-36s: %s", leaf->GetName(), descr.data());
         else {
            // remove leaf name from the description
            std::string title(leaf->GetTitle());
            size_t ind = title.find('[');
            if (ind >= title.npos)
               FATAL("string::find('[') failed");

            Printf("  %-36s: %s%s", leaf->GetName(), descr.data(), &title[ind]);
         }
      }

   } // leaf loop
}

//______________________________________________________________________________
void TreeReader::GetEntry(Long64_t entry)
{
   /* Sets event number to retrieve next time TreeReader::Get*() called.
    */

   if (fTree->IsA() != TChain::Class())
      fEntry = entry;

   // TChain requires special treatment
   else {
      fEntry = ((TChain*)fTree)->LoadTree(entry);

      // reset caches on switching to new TTree
      if (fTreeNum != ((TChain*)fTree)->GetTreeNumber()) {
         fLeafIdx.clear();
         fLeafType.clear();
         fLeafAddr.clear();
         fLeafValue.clear();

         fTreeNum = ((TChain*)fTree)->GetTreeNumber();
      }
   }

   // reset cache of addresses to leaf payloads
   for (size_t i = 0; i < fLeafValue.size(); i++)
      fLeafValue[i] = NULL;
}

//______________________________________________________________________________
void* TreeReader::GetPtr(const char* branch_name, ETypes cktype, Int_t* nsize)
{
   /* Returns either a pointer to first element of one- or multi-dimensional
    * array or a pointer to an unprocessed object.
    *
    * In particular, for vector<...> arrays of elementary data types as well as
    * for vector<vector<int> > and vector<vector<float> > arrays,
    * vector::front() is returned. In this way, the code will work both with old
    * ntuples as well as with the new ntuples (provided that underlying data
    * types are the same in both cases).
    *
    * branch_name = name of branch to search for;
    * cktype = if not kVoidPtr, an additional type verification is performed;
    * nsize = if not NULL, filled with vector::size() for vector<...> branches.
    */

   // entry number in fLeafAddr, fLeafType and fLeafValue
   int i;

   // find the entry number
   std::map<std::string,int>::const_iterator got = fLeafIdx.find(branch_name);

   if (got != fLeafIdx.end())
      i = got->second;
   else {
      // this code is executed once per branch and per root file
      FindLeaf(branch_name);
      i = fLeafType.size() - 1;
   }

   // verify leaf type, if requested
   if (cktype != kVoidPtr) {
      if (cktype == kArrFloat) {
         if (fLeafType[i] != kArrFloatTLeaf && fLeafType[i] != kArrFloatVector)
            FATAL(Form("branch is not of type Float_t*: %s", branch_name));
      } else if (cktype == kArrInt) {
         if (fLeafType[i] != kArrIntTLeaf && fLeafType[i] != kArrIntVector &&
             fLeafType[i] != kArrUIntVector)
            FATAL(Form("branch is not of type (U)Int_t*: %s", branch_name));
      } else if (cktype == kArrChar) {
         if (fLeafType[i] != kArrCharTLeaf && fLeafType[i] != kArrCharVector &&
             fLeafType[i] != kArrUCharVector)
            FATAL(Form("branch is not of type (U)Char_t*: %s", branch_name));
      } else if (cktype == kArrShort) {
         if (fLeafType[i] != kArrShortTLeaf && fLeafType[i] != kArrShortVector &&
             fLeafType[i] != kArrUShortVector)
            FATAL(Form("branch is not of type (U)Short_t*: %s", branch_name));
      } else if (cktype == kArrLong64) {
         if (fLeafType[i] != kArrLong64TLeaf && fLeafType[i] != kArrLong64Vector &&
             fLeafType[i] != kArrULong64Vector)
            FATAL(Form("branch is not of type (U)Long64_t*: %s", branch_name));
      } else if (cktype == kTObject) {
         if (fLeafType[i] != kTObject)
            FATAL(Form("branch content is not inherited from TObject: %s", branch_name));
      } 
      else if (cktype == kArrString) {
	if (fLeafType[i] != kArrStringVector)
	   FATAL(Form("branch is not of type string*: %s", branch_name));
      }      
      else
         if (cktype != fLeafType[i])
            FATAL(Form("invalid branch type requested: %s", branch_name));
   }

   // load branch contents into memory, if necessary
   if (!fLeafValue[i]) {
      // NOTE: it is assumed here that the corresponding branch is not disabled
      TBranch* br = fLeafAddr[i]->GetBranch();
      if (!br)
         FATAL(Form("TLeaf::GetBranch() failed: %s", branch_name));
      if (br->GetEntry(fEntry) < 0)
         FATAL(Form("TBranch::GetEntry() failed: %s", branch_name));

      // pointer to actual leaf payload
      void* ptr = fLeafAddr[i]->GetValuePointer();
      if (fLeafType[i] == kTObject)
         ptr = *((void**)ptr);

      // cache address to payload
      if (fLeafType[i] == kArrFloatVector)
         fLeafValue[i] = &((std::vector<float>*)ptr)->front();
      else if (fLeafType[i] == kArrIntVector)
         fLeafValue[i] = &((std::vector<int>*)ptr)->front();
      else if (fLeafType[i] == kArrUIntVector)
         fLeafValue[i] = &((std::vector<unsigned int>*)ptr)->front();
      else if (fLeafType[i] == kArrCharVector)
         fLeafValue[i] = &((std::vector<char>*)ptr)->front();
      else if (fLeafType[i] == kArrUCharVector)
         fLeafValue[i] = &((std::vector<unsigned char>*)ptr)->front();
      else if (fLeafType[i] == kArrShortVector)
         fLeafValue[i] = &((std::vector<short>*)ptr)->front();
      else if (fLeafType[i] == kArrUShortVector)
         fLeafValue[i] = &((std::vector<unsigned short>*)ptr)->front();
      else if (fLeafType[i] == kArrLong64Vector)
         fLeafValue[i] = &((std::vector<long>*)ptr)->front();
      else if (fLeafType[i] == kArrULong64Vector)
         fLeafValue[i] = &((std::vector<unsigned long>*)ptr)->front();
      else if (fLeafType[i] == kArrStringVector)
	{
	  fLeafValue[i] = &((std::vector<std::string>*)ptr)->front();
	  fStringVectorSize = ((std::vector<std::string>*)ptr)->size();
	}
      else
         fLeafValue[i] = ptr;
   }

   // special case of vector<vector<float> > and vector<vector<int> >
   if (fLeafType[i] == kArrVectFloat) {
      if (nsize)
         *nsize = (Int_t) ((std::vector<std::vector<float> >*)fLeafValue[i])->size();
      return &((std::vector<std::vector<float> >*)fLeafValue[i])->front();
   }
   else if (fLeafType[i] == kArrVectInt) {
      if (nsize)
         *nsize = (Int_t) ((std::vector<std::vector<int> >*)fLeafValue[i])->size();
      return &((std::vector<std::vector<int> >*)fLeafValue[i])->front();
   }

   return fLeafValue[i];
}

//______________________________________________________________________________
void TreeReader::InitSingleTTree(const char* path, const char* treename)
{
   /* Common code for the class constructors: open single TTree.
    */

   // verify sizes of elementary data types
   if (sizeof(int) != 4 || sizeof(long) != 8 || sizeof(float) != 4)
      FATAL("int/long/float of unsupported size");

   // current working directory
   TDirectory* wd = gDirectory;

   // open file with TTree
   fFile = TFile::Open(path);
   if (!fFile || fFile->IsZombie())
      FATAL("TFile::Open() failed");

   // cd back into previous current working directory
   if (wd) wd->cd();
   else gDirectory = 0;

   // get tree
   fTree = dynamic_cast<TTree*>(fFile->Get(treename));
   if (!fTree)
      FATAL(Form("TTree \"%s\" not found:", treename));

   // be 100% sure: check explicitly object's class
   if (((TObject*)fTree)->IsA() != TTree::Class())
      FATAL(Form("\"%s\" is not a TTree", treename));

   // find out availability of MC truth info (check existence of "nMC" branch)
   fkMC = fTree->GetBranch("nMC") ? kTRUE : kFALSE;
}

//______________________________________________________________________________
void TreeReader::InitTChain(const char** paths, int npaths, const char* treename)
{
   /* Common code for the class constructors: make TChain of several TTrees.
    */

   // verify sizes of elementary data types
   if (sizeof(short) != 2 || sizeof(int) != 4 || sizeof(long) != 8 || sizeof(float) != 4)
      FATAL("short/int/long/float of unsupported size");

   // be smart
   if (npaths == 1) {
      InitSingleTTree(paths[0], treename);
      return;
   }

   fTree = new TChain(treename);

   // add root files with TTrees, reading the number of entries in each file
   for (int i = 0; i < npaths; i++)
      if (((TChain*)fTree)->AddFile(paths[i], 0) != 1)
         FATAL("TChain::AddFile() failed");

   // find out availability of MC truth info (check existence of "nMC" branch)
   fkMC = fTree->GetBranch("nMC") ? kTRUE : kFALSE;
}

//______________________________________________________________________________
void TreeReader::FindLeaf(const char* bname)
{
   /* Finds leaf, determines type of its contents and fills the cache of leafs
    * accordingly.
    */

   // find leaf
   TLeaf* leaf = fTree->FindLeaf(bname);
   if (!leaf)
      FATAL(Form("leaf not found: %s", bname));

   // actual data type of the leaf's payload
   ETypes type;

   // vector<...> arrays and other objects not inherited from TObject
   if (leaf->IsA() == TLeafElement::Class()) {
      std::string descr(leaf->GetBranch()->GetClassName());

      if (descr.compare("vector<float>") == 0 ||
          descr.compare("vector<Float_t>") == 0)
         type = kArrFloatVector;
      else if (descr.compare("vector<int>") == 0 ||
               descr.compare("vector<signed>") == 0 ||
               descr.compare("vector<signed int>") == 0 ||
               descr.compare("vector<Int_t>") == 0)
         type = kArrIntVector;
      else if (descr.compare("vector<unsigned>") == 0 ||
               descr.compare("vector<unsigned int>") == 0 ||
               descr.compare("vector<UInt_t>") == 0)
         type = kArrUIntVector;
      else if (descr.compare("vector<char>") == 0 ||
               descr.compare("vector<Char_t>") == 0)
         type = kArrCharVector;
      else if (descr.compare("vector<unsigned char>") == 0 ||
               descr.compare("vector<UChar_t>") == 0)
         type = kArrUCharVector;
      else if (descr.compare("vector<short>") == 0 ||
               descr.compare("vector<short int>") == 0 ||
               descr.compare("vector<signed short>") == 0 ||
               descr.compare("vector<signed short int>") == 0 ||
               descr.compare("vector<Short_t>") == 0)
         type = kArrShortVector;
      else if (descr.compare("vector<unsigned short>") == 0 ||
               descr.compare("vector<unsigned short int>") == 0 ||
               descr.compare("vector<UShort_t>") == 0)
         type = kArrUShortVector;
      else if (descr.compare("vector<long>") == 0 ||
               descr.compare("vector<long int>") == 0 ||
               descr.compare("vector<signed long>") == 0 ||
               descr.compare("vector<signed long int>") == 0 ||
               descr.compare("vector<Long64_t>") == 0 ||
               descr.compare("vector<Long>") == 0)
         type = kArrLong64Vector;
      else if (descr.compare("vector<unsigned long>") == 0 ||
               descr.compare("vector<unsigned long int>") == 0 ||
               descr.compare("vector<ULong64_t>") == 0 ||
               descr.compare("vector<ULong_t>") == 0)
         type = kArrULong64Vector;
      else if (descr.compare("vector<vector<float> >") == 0 ||
               descr.compare("vector<vector<Float_t> >") == 0)
         type = kArrVectFloat;
      else if (descr.compare("vector<vector<int> >") == 0 ||
               descr.compare("vector<vector<signed> >") == 0 ||
               descr.compare("vector<vector<signed int> >") == 0 ||
               descr.compare("vector<vector<Int_t> >") == 0)
         type = kArrVectInt;
     else if (descr.compare("vector<string>") == 0)
         type = kArrStringVector;
     else
         type = kVoidPtr;
   } // TLeafElement

   // objects inherited from TObject
   else if (leaf->IsA() == TLeafObject::Class())
      type = kTObject;

   // fixed/variable length arrays of elementary data types
   else if (leaf->GetLeafCount() || leaf->GetLenStatic() > 1) {
      if (leaf->IsA() == TLeafF::Class())
         type = kArrFloatTLeaf;
      else if (leaf->IsA() == TLeafI::Class())
         type = kArrIntTLeaf;
      else if (leaf->IsA() == TLeafB::Class())
         type = kArrCharTLeaf;
      else if (leaf->IsA() == TLeafS::Class())
         type = kArrShortTLeaf;
      else if (leaf->IsA() == TLeafL::Class())
         type = kArrLong64TLeaf;
      else if (leaf->IsA() == TLeafO::Class())
         type = kArrBool;
      else
         FATAL(Form("branch contains an unknown data type: %s", bname));
   }

   // single variables of elementary data types
   else if (!leaf->GetLeafCount() && leaf->GetLenStatic() == 1) {
      if (leaf->IsA() == TLeafF::Class())
         type = kFloat;
      else if (leaf->IsA() == TLeafI::Class())
         type = kInt;
      else if (leaf->IsA() == TLeafB::Class())
         type = kChar;
      else if (leaf->IsA() == TLeafS::Class())
         type = kShort;
      else if (leaf->IsA() == TLeafD::Class())
         type = kDouble;
      else if (leaf->IsA() == TLeafL::Class())
         type = kLong64;
      else if (leaf->IsA() == TLeafO::Class())
         type = kBool;
      else
         FATAL(Form("branch contains an unknown data type: %s", bname));
   }

   // unknown TLeaf variant
   else
      FATAL(Form("branch contains an unknown data type: %s", bname));

   // update the cache
   fLeafIdx[bname] = (int)fLeafAddr.size();
   fLeafAddr.push_back(leaf);
   fLeafType.push_back(type);
   fLeafValue.push_back(NULL);
}

#endif
