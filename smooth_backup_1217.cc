/* ========================================================================= *
 *                                                                           *
 *                               OpenMesh                                    *
 *           Copyright (c) 2001-2015, RWTH-Aachen University                 *
 *           Department of Computer Graphics and Multimedia                  *
 *                          All rights reserved.                             *
 *                            www.openmesh.org                               *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * This file is part of OpenMesh.                                            *
 *---------------------------------------------------------------------------*
 *                                                                           *
 * Redistribution and use in source and binary forms, with or without        *
 * modification, are permitted provided that the following conditions        *
 * are met:                                                                  *
 *                                                                           *
 * 1. Redistributions of source code must retain the above copyright notice, *
 *    this list of conditions and the following disclaimer.                  *
 *                                                                           *
 * 2. Redistributions in binary form must reproduce the above copyright      *
 *    notice, this list of conditions and the following disclaimer in the    *
 *    documentation and/or other materials provided with the distribution.   *
 *                                                                           *
 * 3. Neither the name of the copyright holder nor the names of its          *
 *    contributors may be used to endorse or promote products derived from   *
 *    this software without specific prior written permission.               *
 *                                                                           *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS       *
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED *
 * TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A           *
 * PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER *
 * OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,  *
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,       *
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR        *
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF    *
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING      *
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS        *
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.              *
 *                                                                           *
 * ========================================================================= */

#include <iostream>
#define no_init_all
#include <Windows.h>
 // -------------------- OpenMesh
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>
 // -------------------- OpenMesh

 //--------------------- allquadrics
#include "quadricfitting.h"
#include "view.h"
#include <fstream>
#include <map>
 //--------------------- allquadrics
#include <numeric>

using namespace allquadrics;
using namespace std;

extern "C" { FILE __iob_func[3] = { *stdin,*stdout,*stderr }; }

// ----------------------------------------------------------------------------

typedef OpenMesh::PolyMesh_ArrayKernelT<>  MyMesh;



// ----------------------------------------------------------------------------
// Build a simple cube and write it to std::cout

void outputQuadric(allquadrics::Quadric &q) {
	cout << q.q[0];
	for (int ii = 1; ii < 10; ii++) {
		cout << ", " << q.q[ii];
	}
}

// ----------------------------------------------------------------------------
// 사이각을 구하기 위한 함수
#define M_PI 3.14159

double angle(vec3 a, vec3 b)
{
	double rad = acos(a * b / (a.length() * b.length()));

	double deg = rad * 180 / M_PI;

	return deg;
}
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------


bool sorting_points(MyMesh mesh, MyMesh::Point a, MyMesh::Point b)
{
	return true;
}

// ----------------------------------------------------------------------------


int main(int argc, char* argv[])
{

	//==============================================openmesh==============================================//
	MyMesh mesh;

	// generate vertices

	MyMesh::VertexHandle vhandle[9];

	vhandle[0] = mesh.add_vertex(MyMesh::Point(-1, -1, 1));
	vhandle[1] = mesh.add_vertex(MyMesh::Point(1, -1, 1));
	vhandle[2] = mesh.add_vertex(MyMesh::Point(1, 1, 1));
	vhandle[3] = mesh.add_vertex(MyMesh::Point(-1, 1, 1));
	vhandle[4] = mesh.add_vertex(MyMesh::Point(-1, -1, -1));
	vhandle[5] = mesh.add_vertex(MyMesh::Point(1, -1, -1));
	vhandle[6] = mesh.add_vertex(MyMesh::Point(1, 1, -1));
	vhandle[7] = mesh.add_vertex(MyMesh::Point(-1, 1, -1));

	// generate (quadrilateral) faces

	std::vector<MyMesh::VertexHandle>  face_vhandles;

	face_vhandles.clear();
	face_vhandles.push_back(vhandle[0]);
	face_vhandles.push_back(vhandle[1]);
	face_vhandles.push_back(vhandle[2]);
	face_vhandles.push_back(vhandle[3]);
	mesh.add_face(face_vhandles);


	face_vhandles.clear();
	face_vhandles.push_back(vhandle[7]);
	face_vhandles.push_back(vhandle[6]);
	face_vhandles.push_back(vhandle[5]);
	face_vhandles.push_back(vhandle[4]);
	mesh.add_face(face_vhandles);

	face_vhandles.clear();
	face_vhandles.push_back(vhandle[1]);
	face_vhandles.push_back(vhandle[0]);
	face_vhandles.push_back(vhandle[4]);
	face_vhandles.push_back(vhandle[5]);
	mesh.add_face(face_vhandles);

	face_vhandles.clear();
	face_vhandles.push_back(vhandle[2]);
	face_vhandles.push_back(vhandle[1]);
	face_vhandles.push_back(vhandle[5]);
	face_vhandles.push_back(vhandle[6]);
	mesh.add_face(face_vhandles);

	face_vhandles.clear();
	face_vhandles.push_back(vhandle[3]);
	face_vhandles.push_back(vhandle[2]);
	face_vhandles.push_back(vhandle[6]);
	face_vhandles.push_back(vhandle[7]);
	mesh.add_face(face_vhandles);

	/*
	face_vhandles.clear();
	face_vhandles.push_back(vhandle[0]);
	face_vhandles.push_back(vhandle[3]);
	face_vhandles.push_back(vhandle[7]);
	face_vhandles.push_back(vhandle[4]);
	mesh.add_face(face_vhandles);
	*/

	//==============================================openmesh==============================================//

	//===============================================Save=================================================//
	// write mesh to output.obj
	try
	{
		if (!OpenMesh::IO::write_mesh(mesh, "output2.obj"))
		{
			std::cerr << "Cannot write mesh to file 'output.off'" << std::endl;
			return 1;
		}
	}
	catch (std::exception& x)
	{
		std::cerr << x.what() << std::endl;
		return 1;
	}

	//===============================================Save=================================================//

		// boundary 찾기
	Sleep(1000);

	MyMesh mesh_findBoundary;
	if (!OpenMesh::IO::read_mesh(mesh_findBoundary, "../Smoothing/output2.obj"))
	{
		std::cerr << "Error: Cannot read mesh from " << std::endl;
		return 1;
	}
	else
	{
		std::cout << "read mesh complete" << std::endl;
	}

	vector<MyMesh::VertexHandle> vhandle_boundary;
	vector<MyMesh::FaceHandle> vhandle_boundary_f;
	vector<MyMesh::HalfedgeHandle> vhandle_boundary_h;
	vector<MyMesh::EdgeHandle> vhandle_boundary_e;
	OpenMesh::PolyConnectivity PC;

	MyMesh::VertexIter          v_it, v_end(mesh_findBoundary.vertices_end());
	MyMesh::VertexVertexIter    vv_it;
	MyMesh::Point               cog;


	MyMesh::FaceIter	f_it, f_end(mesh_findBoundary.faces_end());
	MyMesh::HalfedgeIter	h_it, h_end(mesh_findBoundary.halfedges_end());
	MyMesh::EdgeIter	e_it, e_end(mesh_findBoundary.edges_end());

	std::vector<data_pnw> data;

	int v, f, h, e;

	v = f = h = e = 0;

	for (v_it = mesh_findBoundary.vertices_begin(); v_it != v_end; ++v_it)
	{
		if (mesh_findBoundary.is_boundary(v_it))
		{
			cout << "vertex 하나 넣었음" << endl;

			MyMesh::Point p = mesh_findBoundary.point(v_it);

			cout << p[0] << "\t" << p[1] << "\t" << p[2] << endl;

			vhandle_boundary.push_back(v_it);

			if (!vhandle_boundary.empty())
			{
				MyMesh::HalfedgeHandle hh = mesh_findBoundary.halfedge_handle(v_it);
				MyMesh::VertexHandle vh = mesh_findBoundary.to_vertex_handle(hh);
				MyMesh::HalfedgeHandle hh2 = mesh_findBoundary.halfedge_handle(vh);
				MyMesh::VertexHandle vh2 = mesh_findBoundary.to_vertex_handle(hh2);
				cout << "hh.idx : " << hh.idx() << endl;
				cout << "vh.idx : " << vh.idx() << endl;
				cout << "hh2.idx : " << hh2.idx() << endl;
				cout << "vh2.idx : " << vh2.idx() << endl;
				cout << endl;
			}

			v++;
		}
	}
	/*
	for (f_it = mesh_findBoundary.faces_begin(); f_it != f_end; ++f_it)
	{
		if (mesh_findBoundary.is_boundary(f_it))
		{
			cout << "face 하나 넣었음" << endl;
			vhandle_boundary_f.push_back(f_it);
			f++;
		}
	}
	for (h_it = mesh_findBoundary.halfedges_begin(); h_it != h_end; ++h_it)
	{
		if (mesh_findBoundary.is_boundary(h_it))
		{
			cout << "halfedge 하나 넣었음" << endl;
			vhandle_boundary_h.push_back(h_it);
			h++;
		}
	}
	for (e_it = mesh_findBoundary.edges_begin(); e_it != e_end; ++e_it)
	{
		if (mesh_findBoundary.is_boundary(e_it))
		{
			cout << "edge 하나 넣었음" << endl;
			vhandle_boundary_e.push_back(e_it);
			e++;
		}
	}

	cout << "vertex : " << v << "개" << endl;
	cout << "face : " << f << "개" << endl;
	cout << "halfedge : " << h << "개" << endl;
	cout << "edge : " << e << "개" << endl;
	*/
	
	MyMesh outputMesh;

	MyMesh::VertexHandle vhandle_output[4];

	// vertex handle로 채우기

	vhandle_output[0] = outputMesh.add_vertex(mesh_findBoundary.point(vhandle_boundary[0]));
	vhandle_output[1] = outputMesh.add_vertex(mesh_findBoundary.point(vhandle_boundary[1]));
	vhandle_output[2] = outputMesh.add_vertex(mesh_findBoundary.point(vhandle_boundary[2]));
	vhandle_output[3] = outputMesh.add_vertex(mesh_findBoundary.point(vhandle_boundary[3]));

	face_vhandles.clear();
	face_vhandles.push_back(vhandle_output[2]);
	face_vhandles.push_back(vhandle_output[0]);
	face_vhandles.push_back(vhandle_output[1]);
	face_vhandles.push_back(vhandle_output[3]);
	outputMesh.add_face(face_vhandles);


	vector<MyMesh::Point> holes;

	/*
	for (int i = 0; i < vhandle_boundary.size(); i++)
	{
		holes.push_back(mesh_findBoundary.point(vhandle_boundary[i]));
	}
	*/

	holes.push_back(mesh_findBoundary.point(vhandle_boundary[0]));
	auto vertex_startPoint = vhandle_boundary[0];
	auto vertex_tmp = vhandle_boundary[0];
	MyMesh::HalfedgeHandle hh = mesh_findBoundary.halfedge_handle(vertex_tmp);
	MyMesh::VertexHandle vh = mesh_findBoundary.to_vertex_handle(hh);

	vhandle_boundary.erase(vhandle_boundary.begin());
	
	while (!vhandle_boundary.empty())
	{
		if (vhandle_boundary.size() == 1 && vh.idx() == vertex_startPoint.idx())
		{
			holes.push_back(mesh_findBoundary.point(vertex_tmp));
			break;
		}
		else
		{
			auto p = mesh_findBoundary.point(vertex_tmp);

			cout << "Error : mesh에 포함되지 않은 Vertex가 있습니다. Vertex idx : " << vertex_tmp.idx() << endl;
			cout << "Vertex Coordinate : " << p[0] << " " << p[1] << " " << p[2] << endl;
			return 1;
		}

		for (auto it = vhandle_boundary.begin(); it != vhandle_boundary.end(); it++)
		{
			if (vh.idx == it->idx())
			{
				holes.push_back(mesh_findBoundary.point(*it));
				vertex_tmp = *it;
				
				hh = mesh_findBoundary.halfedge_handle(vertex_tmp);
				vh = mesh_findBoundary.to_vertex_handle(hh);

				vhandle_boundary.erase(it);
				break;
			}
		}		   
		cout << "point로 변환" << endl;
	}

	for (auto i : holes)
	{
		cout << i[0] << "\t" << i[1] << "\t" << i[2] << endl;
	}
	
	vector<double> angles;

	for (int i = 0; i < holes.size(); i++)
	{
		vec3 tmp(holes[i][0], holes[i][1], holes[i][2]);
	
		if (i == 0)
		{
			vec3 a(holes[holes.size() - 1][0], holes[holes.size() - 1][1], holes[holes.size() - 1][2]);
			vec3 b(holes[i + 1][0], holes[i + 1][1], holes[i + 1][2]);

			a = a - tmp;
			b = b - tmp;
			
			angles.push_back(angle(a, b));
		}
		else if (i == holes.size() - 1)
		{
			vec3 a(holes[i - 1][0], holes[i - 1][1], holes[i - 1][2]);
			vec3 b(holes[0][0], holes[0][1], holes[0][2]);

			a = a - tmp;
			b = b - tmp;

			angles.push_back(angle(a, b));
		}
		else
		{
			vec3 a(holes[i - 1][0], holes[i - 1][1], holes[i - 1][2]);
			vec3 b(holes[i + 1][0], holes[i + 1][1], holes[i + 1][2]);

			a = a - tmp;
			b = b - tmp;

			angles.push_back(angle(a, b));
		}
	}

	for (auto i : angles)
	{
		cout << "사이각 : " << i << endl;
	}


	try
	{
		if (!OpenMesh::IO::write_mesh(outputMesh, "output2_cover.obj"))
		{
			std::cerr << "Cannot write mesh to file 'output.off'" << std::endl;
			return 1;
		}
	}
	catch (std::exception& x)
	{
		std::cerr << x.what() << std::endl;
		return 1;
	}
	

	//============================================allquadrics=============================================//
	/*
	vector<allquadrics::Quadric> qfits;
	fitAllQuadricTypes(data, qfits);
	vector<allquadrics::TriangleMesh> meshes; meshes.resize(qfits.size());

	for (int i = 0; i < qfits.size(); i++)
	{
		outputQuadric(qfits[i]);
	}

	cout << endl << endl;
	*/

	//============================================allquadrics=============================================//
}