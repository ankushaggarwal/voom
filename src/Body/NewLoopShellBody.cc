
//----------------------------------------------------------------------
//
//                    William S. Klug, Feng Feng
//                University of California Los Angeles
//                 (C) 2004-2007 All Rights Reserved
//
//----------------------------------------------------------------------

/*! 
  \file NewLoopShellBody.cc

  \brief Implementation of the NewLoopShellBody class, corresponding to
  the concept of a thin shell body.  The body is composed of
  subdivision shell elements

*/

#include <string>
#include <fstream>
#include <ctime>
#include <blitz/array-impl.h>
#include "HalfEdgeMesh.h"
#include "LoopGhostBC.h"

#if defined(_OPENMP)
#include <omp.h>
#endif

namespace voom
{
  using std::endl;
  using std::cout;

  template<class Material_t>
  void NewLoopShellBody<Material_t>::initializeBody
  (Material_t material,
   ConnectivityContainer connectivities,
   const NodeContainer & nodes,
   const unsigned quadOrder) {

    _volume = 0.0;
    _area = 0.0;
  
    // initialize _nDOF and find shell (position) nodes
    _dof = 0;
    _nodes = nodes;
    for(ConstNodeIterator n=_nodes.begin(); n!=_nodes.end(); n++) {
      _dof+=(*n)->dof();
      FeNode_t * sn = dynamic_cast<FeNode_t*>(*n);
      if( sn ) _shellNodes.push_back(sn); 
    }
        
    std::cout << "NewLoopShellBody::initializeBody(): dof = "<< _dof 
	      << std::endl;

    //------------------------------------------------------------------------
    // build elements using simple data structure
    std::cout << "Building elements..." << std::endl;
    clock_t t1=clock();

    HalfEdgeMesh * mesh = new HalfEdgeMesh(connectivities, _shellNodes.size());
    
    // add a boundary layer of ghost faces and vertices for each boundary loop

    std::cout << "Adding ghost vertices and faces along " 
	      << mesh->boundaryLoops.size() << " boundary loops..." 
	      << std::endl;
    int numberGhostVertices=0;
    int numberGhostFaces=0;
    for( int L=0; L<mesh->boundaryLoops.size(); L++ ) {
//       std::cout << "Loop " << L << " of " << mesh->boundaryLoops.size() 
// 		<< std::endl;

      std::vector<HalfEdge*> & boundaryEdges = mesh->boundaryLoops[L];

      // step 1: add a ghost node and face for each boundary edge,
      // also create ghostBC for each boundary edge
      for(int h=0; h<boundaryEdges.size(); h++) {

	HalfEdge * H = boundaryEdges[h];
	
// 	std::cout << "\tH=" << H->id ;

	// add new ghost node
	int V = _shellNodes.size();
	NodeBase::DofIndexMap idx(3,-1);
	Vector3D X(0.0);
	FeNode_t * N = new FeNode_t(V, idx, X, X);
	_shellNodes.push_back(N);
	numberGhostVertices++;
	
	// Create ghost triangle V0 V2 V and ghostBC for adjacent faces
	//      
	//          V
	//         / \
	//       /     \
	//     V0-------V2
	//       \  H  /
	// H->next \ /  H->prev
	//          V1
	//
	int V0 = H->vertex->id;
	int V1 = H->next->vertex->id;
	int V2 = H->prev->vertex->id;

	FeNode_t * N0 = _shellNodes[ V0 ];
	FeNode_t * N1 = _shellNodes[ V1 ];
	FeNode_t * N2 = _shellNodes[ V2 ];

	ElementConnectivity c;
	c = V0, V2, V;
	connectivities.push_back(c);
	numberGhostFaces++;
	
// 	std::cout << "\tF=" << connectivities.size()-1 
// 		  << "\t(" << V0 
// 		  << ", "  << V2
// 		  << ", "  << V
// 		  << ")" << std::endl;


	LoopGhostBC * bc = new LoopGhostBC(N0,N1,N2,N);
	_constraints.push_back( bc );
      
      }

      // step 2: add ghost faces
      // connecting neighboring ghost vertices
      for(int h=0; h<boundaryEdges.size(); h++) {
	//            
	//      VHH-------VH
	//       / \ Add / \  
	//     /     \ /     \
	//    *-------V-------*
	//        HH      H
	//
	HalfEdge * H = boundaryEdges[h];
	HalfEdge * HH = boundaryEdges[(h+1) % boundaryEdges.size()];

	int V = H->vertex->id;
	int VH = _shellNodes.size() - boundaryEdges.size() + h;
	int VHH = _shellNodes.size() - boundaryEdges.size() + (h+1) % boundaryEdges.size();
// 	int VH = mesh->vertices.size() + h;
// 	int VHH = mesh->vertices.size() + (h+1) % boundaryEdges.size();

	ElementConnectivity c;
	c = V, VH, VHH;
	connectivities.push_back(c);
	numberGhostFaces++;
// 	std::cout << "\t\tF=" << connectivities.size()-1 
// 		  << "\t(" << V 
// 		  << ", "  << VH
// 		  << ", "  << VHH
// 		  << ")" << std::endl;
      }
    }

    std::cout << "Added " << numberGhostVertices << " ghost vertices and " 
	      << numberGhostFaces << " ghost faces." << std::endl;

    std::cout << "Rebuilding halfedge mesh"<<std::endl;
    std::cout.flush();

    // rebuild halfedge mesh
    delete mesh;
    mesh = new HalfEdgeMesh(connectivities, _shellNodes.size());
    
    std::cout << "." << std::endl;

//     // step 1: assuming a single boundary, walk around it starting
//     // from one edge

//     // find a boundary edge to start
//     HalfEdge * Hstart=0;
//     for(int h=0; h<mesh->halfEdges.size(); h++) {
//       HalfEdge * H = mesh->halfEdges[h];
//       if( H->opposite == 0 ) {
// 	Hstart = H; 
// 	break;
//       }
//     }

//     // found one, now walk around boundary, storing boundary edges in
//     // order
//     std::vector< HalfEdge * > boundaryEdges;
//     if( Hstart !=0 ) {
//       HalfEdge * H = Hstart; 
//       do {
// 	boundaryEdges.push_back( H );
// 	// find next boundary edge (CCW around the boundary) by
// 	// walking around H's vertex CW
// 	H = H->next;
// 	while ( H->opposite != 0 ) {
// 	  H = H->opposite->next;
// 	}
//       } while ( H != Hstart );
//     }

//     std::cout << "Identified " << boundaryEdges.size() << " boundary edges." 
// 	      << std::endl;
    
//     if( boundaryEdges.size() > 0 ) {
//       // step 2: add a ghost node and face for each boundary edge,
//       // also create ghostBC for each boundary edge
//       for(int h=0; h<boundaryEdges.size(); h++) {
// 	HalfEdge * H = boundaryEdges[h];
	
// 	// add new ghost node
// 	int V = _shellNodes.size();
// 	NodeBase::DofIndexMap idx(3,-1);
// 	Vector3D X(0.0);
// 	FeNode_t * N = new FeNode_t(V, idx, X, X);
// 	_shellNodes.push_back(N);
	
// 	// Create ghost triangle V0 V2 V and ghostBC for adjacent faces
// 	//      
// 	//          V
// 	//         / \
// 	//       /     \
// 	//     V0-------V2
// 	//       \  H  /
// 	// H->next \ /  H->prev
// 	//          V1
// 	//
// 	int V0 = H->vertex->id;
// 	int V1 = H->next->vertex->id;
// 	int V2 = H->prev->vertex->id;

// 	FeNode_t * N0 = _shellNodes[ V0 ];
// 	FeNode_t * N1 = _shellNodes[ V1 ];
// 	FeNode_t * N2 = _shellNodes[ V2 ];

// 	ElementConnectivity c;
// 	c = V0, V2, V;
// 	connectivities.push_back(c);

// 	LoopGhostBC * bc = new LoopGhostBC(N0,N1,N2,N);
// 	_constraints.push_back( bc );
      
//       }

//       // step 3: add ghost faces
//       // connecting neighboring ghost vertices
//       for(int h=0; h<boundaryEdges.size(); h++) {
// 	//            
// 	//      VHH-------VH
// 	//       / \ Add / \  
// 	//     /     \ /     \
// 	//    *-------V-------*
// 	//        HH      H
// 	//
// 	HalfEdge * H = boundaryEdges[h];
// 	HalfEdge * HH = boundaryEdges[(h+1) % boundaryEdges.size()];

// 	int V = H->vertex->id;
// 	int VH = mesh->vertices.size() + h;
// 	int VHH = mesh->vertices.size() + (h+1) % boundaryEdges.size();

// 	ElementConnectivity c;
// 	c = V, VH, VHH;
// 	connectivities.push_back(c);
//       }

//       // rebuild halfedge mesh
//       delete mesh;
//       mesh = new HalfEdgeMesh(connectivities, _shellNodes.size());
//     }
       
    for(int f=0; f<mesh->faces.size(); f++) {
      Face * F = mesh->faces[f];
      // get valences of corners and find one-ring of
      // next-nearest-neighbor vertices
      typename LoopShell<Material_t>::NodeContainer nds;
      CornerValences v;
      //don't create element for this face if any of the vertices is on the boundary
      if( F->halfEdges[0]->vertex->boundary ) continue;
      if( F->halfEdges[1]->vertex->boundary ) continue;
      if( F->halfEdges[2]->vertex->boundary ) continue;

      for(int a=0; a<3; a++) {
	// add corner vertices
	HalfEdge * H = F->halfEdges[a];
	nds.push_back( _shellNodes[H->vertex->id] );
	// valence is the number of incident halfedges
	v(a) = H->vertex->halfEdges.size();
      }
      // walk around the one-ring adding vertices
      //
      // Start by walking CCW through the vertices incident to the
      // second corner vertex.  Stop at first vertex incident to the
      // third corner vertex.
      int v1=3;
      for(HalfEdge * H = F->halfEdges[1]->opposite->next; 
	  H->vertex->id != F->halfEdges[2]->opposite->next->vertex->id; 
	  H = H->next->opposite->next  ) {
	nds.push_back( _shellNodes[H->vertex->id] );
	v1++;
      }
      // Now continue CCW through the vertices incident to the third
      // corner vertex. Stop at first vertex incident to the first
      // corner vertex.
      int v2=3;
      for(HalfEdge * H = F->halfEdges[2]->opposite->next; 
	  H->vertex->id != F->halfEdges[0]->opposite->next->vertex->id; 
	  H = H->next->opposite->next  ) {
	nds.push_back( _shellNodes[H->vertex->id] );
	v2++;
      }
      // Now continue CCW through the vertices incident to the first
      // corner vertex. Stop at first vertex incident to the second
      // corner vertex.
      int v0=3;
      for(HalfEdge * H = F->halfEdges[0]->opposite->next; 
	  H->vertex->id != F->halfEdges[1]->opposite->next->vertex->id; 
	  H = H->next->opposite->next  ) {
	nds.push_back( _shellNodes[H->vertex->id] );
	v0++;
      }
      assert(v0==v(0));
      assert(v1==v(1));
      assert(v2==v(2));
      const unsigned npn = v(0) + v(1) + v(2) - 6;
      
      FeElement_t * elem 
	= new FeElement_t( TriangleQuadrature(quadOrder), 
			   material,
			   nds,
			   v);
      _shells.push_back(elem);

      _active.push_back(true);
      
    }

    delete mesh;

    _elements.reserve(_shells.size());
    for(ConstFeElementIterator e=_shells.begin(); e!=_shells.end(); e++) {
      _elements.push_back(*e);
    }

    clock_t t2=clock();
    std::cout << "Done building elements.  Elapsed time: "
	      << (double)(t2-t1)/CLOCKS_PER_SEC
	      << std::endl;
    std::cout<< "Number of shells "<<_shells.size()<<"\nNumber of elements "<<_elements.size()<<std::endl;
    // compute mechanics
    compute(false, false, false);
  }

  //! set the spontaneous curvature as the reference configuration curvature
  template< class Material_t >
  void NewLoopShellBody<Material_t>::SetRefMeanCurvature()
  {
    int sBegin=0, sEnd=_shells.size();
    for(int si=sBegin; si<sEnd; si++) {
      if( !_active[si] ) continue;
      FeElement_t* s=_shells[si];
      s->SetRefMeanCurvature();
     }
  }

  //! compute
  template< class Material_t >
  void NewLoopShellBody<Material_t>::compute( bool f0, bool f1, bool f2 )
  { 
    
    int eBegin=0, eEnd=_elements.size();
    int sBegin=0, sEnd=_shells.size();

#ifdef _OPENMP	
#pragma omp parallel for 		\
  schedule(static) default(shared) 
#endif	
    for(int ei=eBegin; ei<eEnd; ei++) {
	if( !_active[ei] ) continue;
	  _elements[ei]->compute( f0, f1, f2 );
    }

    if(f0) { 
// #ifdef _OPENMP	
// #pragma omp parallel for
// #endif	
      for(int ei=eBegin; ei<eEnd; ei++) {
	if( !_active[ei] ) continue;
	Element* e=_elements[ei];
	_energy += e->energy();	
      }
    }
    return;
  }

  template< class Material_t >                                                
  void NewLoopShellBody<Material_t>::cal_curv(std::vector<double> &curv )        
  {                                                                           
    int sBegin=0, sEnd=_shells.size();
    for(int si=sBegin; si<sEnd; si++) {
       if( !_active[si] ) continue;
       FeElement_t* s=_shells[si];
       curv[si]=s->meancurvature();
     }
    return;                                                                   
  }           

  //! create input file used by Paraview, a 3D viewer
  template < class Material_t >
  void NewLoopShellBody< Material_t >::printParaview(const std::string name) const
  {

    int numShells = _shells.size();
    
    blitz::Array<double,1> energy(blitz::shape(numShells));
    blitz::Array<double,1> bendingEnergy(blitz::shape(numShells));
    blitz::Array<double,1> inplaneEnergy(blitz::shape(numShells));
    blitz::Array<double,1> curvature(blitz::shape(numShells));

    energy = 0.0;
    curvature = 0.0;
    bendingEnergy = 0.0;
    inplaneEnergy = 0.0;

    ConstFeElementIterator pe = _shells.begin();
    for ( int e = 0; pe != _shells.end(); pe++, e++) {
      int npts=0;
      typename FeElement_t::ConstQuadPointIterator 
	p = (*pe)->quadraturePoints().begin();
      for( ; p != (*pe)->quadraturePoints().end(); p++){
	energy(e) += p->material.energyDensity();
	bendingEnergy(e) += p->material.bendingEnergy();
	inplaneEnergy(e) += p->material.stretchingEnergy();
	curvature(e) += p->material.meanCurvature();
	npts++;
      }
      energy(e) /= (double)( npts );
      bendingEnergy(e) /= (double)( npts );
      inplaneEnergy(e) /= (double)( npts );
      curvature(e) /= (double)( npts );
    }

    std::string fileName = name + ".vtk";
    std::ofstream ofs(fileName.c_str());
    if (!ofs) {
      std::cout << "can not open output ("
		<< fileName
		<< ") file." << std::endl;
      exit(0);
    }
    
    ////////////////////////////////////////////////////////////////////
    //
    //    Node Section
    //
    ofs << "# vtk DataFile Version 2.0\n"
	<< "Test example" << std::endl
	<< "ASCII" << std::endl
	<< "DATASET POLYDATA" << std::endl
	<< "POINTS  " << _shellNodes.size() << "  double" << std::endl;
    
    //
    // output nodal postions
    ConstFeNodeIterator pn = _shellNodes.begin();
    for ( ; pn!= _shellNodes.end(); pn ++) {
      const Vector3D & nodalPos =  (*pn)->position();
      ofs << std::setprecision(16) 
	  << nodalPos(0) << "  "
	  << nodalPos(1) << "  "
	  << nodalPos(2) << std::endl;
    }

    /////////////////////////////////////////////////////////////////////
    //
    //    Element Section
    //
    int nActiveShells=0;
    for(int e=0; e<_shells.size(); e++) {
      if(_active[e]) nActiveShells++;
    }
    ofs << "POLYGONS  " << nActiveShells << "  "
	<< 4*nActiveShells << std::endl;
    for(int e=0; e<_shells.size(); e++) {
      if(!_active[e])  continue;
      const FeElement_t * pe = _shells[e];
      ofs << 3 << "  "
	  << std::setw(10) << (pe)->nodes()[0] -> id()
	  << std::setw(10) << (pe)->nodes()[1] -> id()
	  << std::setw(10) << (pe)->nodes()[2] -> id()
	  << std::endl;
    }
     	
    //////////////////////////////////////////////////////////////////////////
    //
    //  output color for each element ( corresponding to mean
    //  curvature, energy...)
    //
    ofs << "CELL_DATA    " << nActiveShells << std::endl;
    //
    // output for strain energy
    ofs << "SCALARS    strainEnergy    double    1" << std::endl;
    ofs << "LOOKUP_TABLE default" << std::endl;
    for(int e=0; e<_shells.size(); e++) {
      if(!_active[e])  continue;
      ofs << energy(e) << std::endl;
    }
    ofs << std::endl;

    //
    // output for strain energy
    ofs << "SCALARS    bendingEnergy    double    1" << std::endl;
    ofs << "LOOKUP_TABLE default" << std::endl;
    for(int e=0; e<_shells.size(); e++) {
      if(!_active[e])  continue;
      ofs << bendingEnergy(e) << std::endl;
    }
    ofs << std::endl;

    //
    // output for strain energy
    ofs << "SCALARS    inplaneEnergy    double    1" << std::endl;
    ofs << "LOOKUP_TABLE default" << std::endl;
    for(int e=0; e<_shells.size(); e++) {
      if(!_active[e])  continue;
      ofs << inplaneEnergy(e) << std::endl;
    }
    ofs << std::endl;

    //
    // output for mean curvature
    // since one color is mapping into one element, for
    // several guass points integration, we compute their average value
    ofs << "SCALARS    meanCurvature    double    1" << std::endl;
    ofs << "LOOKUP_TABLE meanCurvature    " << std::endl;		
    for(int e=0; e<_shells.size(); e++) {
      if(!_active[e])  continue;
      ofs << curvature(e) << std::endl;
    }
    ofs << std::endl;

    ofs << endl << "POINT_DATA " << _shellNodes.size() << endl
	<< "VECTORS displacements double" << endl;
    //      << "LOOKUP_TABLE displacements" << endl; 
    //
    // output nodal postions
    pn = _shellNodes.begin();
    for ( ; pn!= _shellNodes.end(); pn ++) {
      Vector3D nodalDisp;
      nodalDisp = (*pn)->point() - (*pn)->position();
      ofs << std::setprecision(16) 
	  << nodalDisp(0)
	  << '\t' <<nodalDisp(1)
	  << '\t' <<nodalDisp(2) << std::endl;
    }
    
    ofs << endl << "VECTORS forces double" << endl;
    //     << "LOOKUP_TABLE displacements" << endl; 
    //
    // output nodal postions
    pn = _shellNodes.begin();
    for ( ; pn!= _shellNodes.end(); pn ++) {
      const Vector3D & nodalForce = (*pn)->force();
      ofs << std::setprecision(16) 
	  << nodalForce(0)
	  << '\t' <<nodalForce(1)
	  << '\t' <<nodalForce(2) << std::endl;
    }

    ofs << endl << "VECTORS contactForces double" << endl;
    //     << "LOOKUP_TABLE displacements" << endl; 
    
    for (int a=0; a < _shellNodes.size(); a++) {
      Vector3D f(0.0);
      for(ConstConstraintIterator c=_constraints.begin(); c!=_constraints.end(); c++) {
	for(int i=0; i<3; i++)	f(i) += (*c)->getForce(a,i);
      }

      ofs << std::setprecision(16) 
	  << f(0) << '\t' <<f(1) << '\t' <<f(2) 
	  << std::endl;
    }
    ofs.close();

    return;
  } 
} // namespace voom
