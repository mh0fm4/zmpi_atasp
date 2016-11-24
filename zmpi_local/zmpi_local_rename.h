/*
 *  Copyright (C) 2012, 2013, 2014, 2015, 2016 Michael Hofmann, Chemnitz University of Technology
 *  
 *  This file is part of the ZMPI All-to-all Specific Library.
 *  
 *  The ZMPI-ATASP is free software: you can redistribute it and/or
 *  modify it under the terms of the GNU Lesser Public License as published
 *  by the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *  
 *  The ZMPI-ATASP is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser Public License for more details.
 *  
 *  You should have received a copy of the GNU Lesser Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */


#ifndef __ZMPI_LOCAL_RENAME_H__
#define __ZMPI_LOCAL_RENAME_H__


#ifndef ZMPI_RENAME

#define ZMPI_RENAME

#define ZMPI_CONCAT(_a_, _b_)           ZMPI_CONCAT_(_a_, _b_)
#define ZMPI_CONCAT_(_a_, _b_)          _a_##_b_

#define ZMPI_CONCONCAT(_a_, _b_, _c_)   ZMPI_CONCONCAT_(_a_, _b_, _c_)
#define ZMPI_CONCONCAT_(_a_, _b_, _c_)  _a_##_b_##_c_

#ifdef ZMPI_PREFIX
# define ZMPI_VAR(_v_)   ZMPI_CONCAT(ZMPI_PREFIX, _v_)
# define ZMPI_FUNC(_f_)  ZMPI_CONCAT(ZMPI_PREFIX, _f_)
#else
# define ZMPI_VAR(_v_)   _v_
# define ZMPI_FUNC(_f_)  _f_
#endif

#endif /* ZMPI_RENAME */


/* zmpil_simple.c */
#define zmpil_simple_create_derived  ZMPI_FUNC(zmpil_simple_create_derived)
#define zmpil_simple_destroy  ZMPI_FUNC(zmpil_simple_destroy)


#endif /* __ZMPI_LOCAL_RENAME_H__ */
