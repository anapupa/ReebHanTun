# "ReebHanTun" software for computing Handle and Tunnel loops of 3D surfaces 

"ReebHanTun" software is developed by the Jyamiti research group headed by 
Prof. Tamal K. Dey at the Department of Computer Science and Engineering 
of The Ohio State University.

The binaries are distributed for: Windows 64bit and Ubuntu Linux 32bit;

## DESCRIPTION


The "ReebHanTun" software can compute two families of non-trivial loops 
called handle and tunnel loops for any orientable 3D surface mesh.
This software is developed based on the following paper:

	* Paper: T. K. Dey, F. Fan and Y. Wang.  
	An Efficient Computation of Handle and Tunnel Loops via Reeb Graphs. 
	ACM SIGGRAPH 2013.
 
"ReebHanTun" software can be used in the following two cases so far:

- Connected, closed and orientable 3D surface meshes:

> For a connected, closed and orientable surface mesh of genus g, 
the "ReebHanTun" software takes the mesh as input and 
computes g Handle loops and g Tunnel loops.
All computed loops are vertex-edge paths on the input mesh.

- Connected and orientable 3D surface meshes with boundaries:

> For a connected and orientable surface mesh of genus g 
with disk type boundaries, 
the "ReebHanTun" software takes the open mesh as input and
computes g Handle loops and g Tunnel loops.
All computed loops are vertex-edge paths on the input non-closed mesh.

---
_**WARNING**_:  
As described in the paper, the "ReebHanTun" software
uses a very simple heuristic to handle meshes with boundaries.
For meshes with small disk type boundaries,  
the "ReebHanTun" software can produce correct output.
However, it is possible that the "ReebHanTun" software will not be 
robust for arbitrary boundaries.  

### USAGE

The usages for two cases are the same.

`./ReebHanTun -I <input_mesh> -O <outfile> [-options]`

```
INPUT PARAMETERS:
<input_mesh>       	: input 3D surface mesh in OFF format:
												OFF
												#v #f #e
												v1_x v1_y v1_z
												v2_x v2_y v2_z
												...
												vn_x vn_y vn_z
												3 v1_1 v1_2 v1_3
												3 v2_1 v2_2 v2_3
												...
												3 vm_1 vm_2 vm_3						

<outfile> : the name that the user gives for the output files

[-options]   : can be one of the following:
	       	-h      (Help information)
	       	-l     	(Display full license) 
```

### OUTPUT FILES:

"ReebHanTun" produces two output files: one mesh file for visualization of 
the Handle / Tunnel loops, and one file that contain the loop information 
of the following format:

(1). "loops_<outfile_prefix>.list" in LIST format 
which can be visualized in the Geomview software
(www.geomview.org).

(2). "loops_<outfile_prefix>.lop": 
the Handle / Tunnel loops file in the following format:

*First, it stores all the Handle loops:

# there are m handle loops 
Handle Loop 0(loop size): a list of vertices on this loop 
Handle Loop 1(loop size): a list of vertices on this loop
	          ...
Handle Loop m-1(loop size): a list of vertices on this loop

*Then, it stores all the Tunnel loops:

# there are m tunnel loops
Tunnel Loop 0(loop size): a list of vertices on this loop
Tunnel Loop 1(loop size): a list of vertices on this loop
	          ...
Tunnel Loop m-1(loop size): a list of vertices on this loop

=================================
EXAMPLE
=================================

Usage example : run "ReebHanTun" for a closed / open "eight.off"  model

Command : "ReebHanTun -I eight.off -O eight"

The input mesh is "eight.off". 
The output files containing handle and tunnel loops are
"loops_eight.list" and "loops_eight.lop"

Command : "ReebHanTun -I eight.off -O eight -h"
The additional '-h' option will display the following 
help information:

ReebHanTun Usage:
  -h                    Help information
  -l                    License information
  -I arg                Input file name
  -O arg                Output file name prefix

Command : "ReebHanTun -I eight.off -O eight -l"
The additional '-l' option will display the license information.

### LEGAL TERMS

THIS SOFTWARE IS PROVIDED "AS-IS". THERE IS NO WARRANTY OF ANY KIND. 
NEITHER THE AUTHORS NOR THE OHIO STATE UNIVERSITY WILL BE LIABLE FOR 
ANY DAMAGES OF ANY KIND, EVEN IF ADVISED OF SUCH POSSIBILITY.

This software was developed (and is copyrighted by) the Jyamiti group at 
The Ohio State University. Please do not redistribute this software. 
This program is for academic research use only. This software uses the 
CGAL library (www.cgal.org), Boost library (www.boost.org) and Ann library
(www.cs.umd.edu/~mount/ANN/) which are covered under their own licenses.

The CGAL library's license 
(which applies to the CGAL library ONLY and NOT to this program itself) is 
as follows:
 
### LICENSE


The CGAL software consists of several parts, each of which is licensed under
an open source license. It is also possible to obtain commercial licenses
from GeometryFactory (www.geometryfactory.com) for all or parts of CGAL.

The source code of the CGAL library can be found in the directories
"src/CGAL", "src/CGALQt" and "include/CGAL". It is specified in each file of
the CGAL library which license applies to it. This is either the GNU Lesser
General Public License (as published by the Free Software Foundation;
version 2.1 of the License) or the Q Public License (version 1.0). The texts
of both licenses can be found in the files LICENSE.LGPL and LICENSE.QPL.  

Distributed along with CGAL (for the users' convenience), but not part of
CGAL, are the following third-party libraries, available under their own
licenses:

- CORE, in the directories "include/CORE" and "src/Core", is licensed under
  the QPL (see LICENSE.QPL).
- OpenNL, in the directory "include/OpenNL", is licensed under the LGPL
  (see include/OpenNL/LICENSE.OPENNL).
- ImageIO, in the directory "examples/Surface_mesher/ImageIO", is licensed
  under the LGPL (see LICENSE.LGPL).

All other files that do not have an explicit copyright notice (e.g., all
examples and some demos) are licensed under a very permissive license. The
exact license text can be found in the file LICENSE.FREE_USE. Note that some
subdirectories have their own copy of LICENSE.FREE_USE. These copies have
the same license text and differ only in the copyright holder.
---------------------------------------------------------------------------

The Boost library's license 
(which applies to the Boost library ONLY and NOT to this program itself) is 
as follows:

LICENSE
---------------------------------------------------------------------------
Boost Software License - Version 1.0 - August 17th, 2003

Permission is hereby granted, free of charge, to any person or organization
obtaining a copy of the software and accompanying documentation covered by
this license (the "Software") to use, reproduce, display, distribute,
execute, and transmit the Software, and to prepare derivative works of the
Software, and to permit third-parties to whom the Software is furnished to
do so, all subject to the following:

The copyright notices in the Software and this entire statement, including
the above license grant, this restriction and the following disclaimer,
must be included in all copies of the Software, in whole or in part, and
all derivative works of the Software, unless such copies or derivative
works are solely in the form of machine-executable object code generated by
a source language processor.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT
SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE
FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE,
ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.
