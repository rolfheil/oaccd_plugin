#
# @BEGIN LICENSE
#
# oaccd by Psi4 Developer, a plugin to:
#
# Psi4: an open-source quantum chemistry software package
#
# Copyright (c) 2007-2017 The Psi4 Developers.
#
# The copyrights for code used from other parties are included in
# the corresponding files.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#
# @END LICENSE
#

cmake_minimum_required(VERSION 3.1 FATAL_ERROR)

project(oaccd CXX)

find_package(psi4 1.1 REQUIRED)

add_psi4_plugin(oaccd oaccd.h oaccd.cc int_trans_rhf.cc mp2_energy.cc ccd_energy.cc biort_libtrans/biort_inttransform.h biort_libtrans/biort_inttransform.cc biort_libtrans/transform_tei.cc biort_libtrans/transform_tei_first_half.cc biort_libtrans/transform_tei_second_half.cc biort_libtrans/biort_transform_moinfo.cc biort_libtrans/biort_fock_construct.cc biort_wfn/biortwfn.h biort_wfn/biortwfn.cc)
