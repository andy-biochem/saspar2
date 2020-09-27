# This file is part of SASPAR2 software package.
# Copyright 2020 Andrei M. Semenov
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#  
#     http://www.apache.org/licenses/LICENSE-2.0
#  
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import numpy


def generate_vectors(q_module, param_n, param_m):
    vectors_ar = []
    dz = 2 / (param_n - 1)
    dfi = (2 * numpy.pi) / param_m

    for i in range(0, param_n):
        z = (dz * i) - 1.0
        teta = numpy.arccos(z)
        for j in range(0, param_m):
            fi = dfi * j
            q_x = q_module * numpy.sin(teta) * numpy.cos(fi)
            q_y = q_module * numpy.sin(teta) * numpy.sin(fi)
            q_z = q_module * numpy.cos(teta)
            vectors_ar.append((q_x, q_y, q_z))

    return vectors_ar