# This file is part of saddle-bags.
#
# saddle-bags is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# saddle-bags is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with saddle-bags. If not, see <http://www.gnu.org/licenses/>.

# This file contains specifications for packaging of saddlebags
# As a standalone executable.  This file is meant to be used with pyinstaller
# http://www.pyinstaller.org/


# -*- mode: python -*-

block_cipher = None

# TODO: I should include the file from punkin chunker.
# It is the GeneRef.fasta file. This is how i sort by gene

a = Analysis(['NanoporeProspectorMain.py'],
             binaries=None,
             datas= [ ('../nit-picker/barcodes/barcodes_96_bc_kit.txt', '') , ('../punkin-chunker/inputData/HLA_ClassI_GeneRef.fasta', '') ],
             hiddenimports=['six', 'packaging', 'packaging.requirements', 'packaging.version', 'packaging.specifiers', 'Tkinter', 'tkFileDialog', 'Tkconstants', 'Bio', 'Bio.Seq', 'Bio.Blast', 'Bio.Blast.NCBIXML', 'matplotlib', 'pylab', 'pysam', 'numpy'],
             hookspath=[],
             runtime_hooks=[],
             excludes=['tkinter'],
             win_no_prefer_redirects=False,
             win_private_assemblies=False,
             cipher=block_cipher)
             
pyz = PYZ(a.pure, a.zipped_data,
             cipher=block_cipher)
             
exe = EXE(pyz,
          a.scripts,
          a.binaries,
          a.zipfiles,
          a.datas,
          name='NanoporeProspectorWindows',
          debug=False,
          strip=False,
          upx=True,
          console=True )

