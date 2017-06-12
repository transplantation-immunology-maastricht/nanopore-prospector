# This file is part of nanopore-prospector.
#
# nanopore-prospector is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# nanopore-prospector is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with nanopore-prospector. If not, see <http://www.gnu.org/licenses/>.

# Version 1.0 

SoftwareVersion = "nanopore-prospector Version 1.0"

import Tkinter
#import sys
from NanoporeProspectorMasterFrame import NanoporeProspectorMasterFrame

if __name__=='__main__':

    # TODO: Read args? What args do i even need?

    root = Tkinter.Tk()
    app = NanoporeProspectorMasterFrame(root)
    root.mainloop()

    print('Done.  Yay.')


