#
#    Copyright (C) 2010-2014 University of Southern California and
#                            Andrew D. Smith
#
#    Authors: Andrew D. Smith
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

all_subdirs=common tools
app_subdirs=tools

all:
	@make -C tools SRC_ROOT=$(CREAD) OPT=1

install:
	@make -C tools SRC_ROOT=$(CREAD) OPT=1 install

clean:
	@for i in $(all_subdirs); do \
		make -C $${i} SRC_ROOT=$(CREAD) clean; \
	done;
