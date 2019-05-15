# encoding: utf-8

# Copyright 2013 Diego Navarro Mellén. All rights reserved.
#
# Redistribution and use in source and binary forms, with or without modification, are
# permitted provided that the following conditions are met:
#
#  1. Redistributions of source code must retain the above copyright notice, this list of
#     conditions and the following disclaimer.
#
#  2. Redistributions in binary form must reproduce the above copyright notice, this list
#     of conditions and the following disclaimer in the documentation and/or other materials
#     provided with the distribution.
#
# THIS SOFTWARE IS PROVIDED BY DIEGO NAVARRO MELLÉN ''AS IS'' AND ANY EXPRESS OR IMPLIED
# WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND
# FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL DIEGO NAVARRO MELLÉN OR
# CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
# ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
# NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
# ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# The views and conclusions contained in the software and documentation are those of the
# authors and should not be interpreted as representing official policies, either expressed
# or implied, of Diego Navarro Mellén.

from __future__ import print_function


# Text attributes
START = '<span '
STYLE = 'style="'
ENDSTYLE = '">'
END   = '</span>'
BOLD = 'font-weight:bold;'
UNDERLINE = 'text-decoration: underline;'
#BLINK = '\033[5m'
#REVERSE = '\033[7m'
#CONCEALED = '\033[7m'

# Foreground colors
FG_BLACK = 'color:black;'
FG_RED = 'color:red;'
FG_GREEN = 'color:green;'
FG_YELLOW = 'color:yellow;'
FG_BLUE = 'color:blue;'
FG_MAGENTA = 'color:magenta;'
FG_CYAN = 'color:cyan;'
FG_WHITE = 'color:white;'

# Background colors
BG_BLACK = 'background-color:black;'
BG_RED = 'background-color:red;'
BG_GREEN = 'background-color:green;'
BG_YELLOW = 'background-color:yellow;'
BG_BLUE = 'background-color:blue;'
BG_MAGENTA = 'background-color:magenta;'
BG_CYAN = 'background-color:cyan;'
BG_WHITE = 'background-color:white;'
