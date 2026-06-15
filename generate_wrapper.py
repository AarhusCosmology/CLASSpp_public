#!/usr/bin/env python
# coding: utf-8

# In[1]:


import argparse
import os
import pathlib

rootdir = pathlib.Path(__file__).parent

_argparser = argparse.ArgumentParser(description='Generate cclassy.pxd from the C++ headers.')
_argparser.add_argument('--output-dir', type=pathlib.Path, default=rootdir,
                        help='Directory to write cclassy.pxd into (default: repo root).')
# parse_known_args: tolerate being run via runpy with foreign argv.
_args, _ = _argparser.parse_known_args()
outdir = _args.output_dir


def read_header_lines(path):
    lines = []
    current_line = ''
    in_block_comment = False
    with open(path, encoding='utf-8') as fid:
        for raw_line in fid:
            raw_line = raw_line.rstrip('\n')
            line = ''
            index = 0
            while index < len(raw_line):
                if in_block_comment:
                    comment_end = raw_line.find('*/', index)
                    if comment_end == -1:
                        index = len(raw_line)
                        continue
                    in_block_comment = False
                    index = comment_end + 2
                    continue
                if raw_line.startswith('//', index):
                    break
                if raw_line.startswith('/*', index):
                    in_block_comment = True
                    index += 2
                    continue
                line += raw_line[index]
                index += 1
            if current_line:
                current_line += ' ' + line.lstrip()
            else:
                current_line = line
            if current_line.endswith('\\'):
                current_line = current_line[:-1]
                continue
            if current_line.strip():
                lines.append(current_line + '\n')
            current_line = ''
    if current_line:
        lines.append(current_line + '\n')
    return lines


def is_comment_only_line(line):
    stripped = line.strip()
    return (not stripped or
            stripped.startswith('//') or
            stripped.startswith('/*') or
            stripped.startswith('/**') or
            stripped.startswith('*') or
            stripped.startswith('*/') or
            stripped.startswith('//@'))

h_files = []
header_dirs = ['include', 'main', 'source', 'species', 'tools']
for d in header_dirs:
    h_files.extend(sorted((rootdir / d).glob('*.h')))


# In[2]:


preample = """
from libcpp cimport bool
from libcpp.map cimport map
from libcpp.memory cimport shared_ptr
from libcpp.pair cimport pair
from libcpp.string cimport string
from libcpp.vector cimport vector

"""


# In[3]:


definition_names = [
    '_MAX_NUMBER_OF_K_FILES_',
    '_Z_PK_NUM_MAX_',
    '_SELECTION_NUM_MAX_',
    '_ERRORMSGSIZE_',
    '_TRUE_',
    '_FALSE_',
    '_SUCCESS_',
    '_FAILURE_',
]
definitions_dict = {}
for file in h_files:
    for line in read_header_lines(file):
        if "#define" not in line:
            continue
        words = line.split()
        if len(words) < 3:
            continue
        if words[0] != "#define":
            continue
        if words[1] in definition_names:
            definitions_dict[words[1]] = words[2]

# In[4]:


enums = []
enums.append('cdef extern from "class.h":')
enums.append('    pair[string, string] get_my_py_error_message()')
enums.append('')
enums.append('    ctypedef char ErrorMsg[_ERRORMSGSIZE_]')
enums.append('')
enum_names = [
    'equation_of_state',
    'file_format',
    'halofit_integral_type',
    'hmcode_baryonic_feedback_model',
    'linear_or_logarithmic',
    'non_linear_method',
    'out_sigmas',
    'pk_outputs',
    'primordial_spectrum_type',
    'selection_type',
    'source_extrapolation',
    'spatial_curvature',
]

for file in h_files:
    enum_found = False
    for line in read_header_lines(file):
        if not enum_found:
            #Check for struct
            if 'enum' not in line:
                continue
            elif '(' in line:
                # Probably a function definition
                continue
            elif '{' not in line:
                # Probably a declaration of an enum variable in a struct.
                continue
            for s in enum_names:
                if 'enum ' + s + ' ' in line or 'enum ' + s + '\n' in line:
                    if '}' not in line:
                        #Assume multi-line enum:
                        enum_found = True
                    line = line.strip()
                    line = line.replace('{',':').replace('}','').replace(';','')
                    enums.append('    ctypedef ' + line)
                    break
        else:
            # Check for multiline enum ended
            if '};' in line:
                enum_found = False
            elif line:
                words = line.strip().split()
                enums.append('       ' + words[0])


# In[5]:


structs = []
struct_names = [
    'background',
    'lensing',
    'nonlinear',
    'output',
    'perturbs',
    'precision',
    'primordial',
    'spectra',
    'thermo',
    'transfers',
]
allowed_types = ['double', 'int', 'short', 'bool', 'std::string']

for file in h_files:
    file_lines = read_header_lines(file)
    struct_found = False
    index_line = 0
    while index_line < len(file_lines):
        line = file_lines[index_line]
        if not struct_found:
            #Check for struct
            if 'struct' not in line:
                index_line += 1
                continue
            for s in struct_names:
                if 'struct ' + s + ' ' in line or 'struct ' + s + '\n' in line:
                    struct_found = True
                    structs.append('    cdef struct ' + s + ':')
                    break
            index_line += 1
        else:
            # Check for struct ended. Match only the struct's closing brace
            # (a line starting with '}'), not any line that merely contains
            # '};' -- e.g. an in-struct array initializer such as
            # `double selection_bias[N] = {1.};` would otherwise be mistaken
            # for the end of the struct and truncate the remaining members.
            if line.strip().startswith('}'):
                structs.append('')
                struct_found = False
                index_line += 1
            else:
                declaration = line.strip()
                while declaration and ';' not in declaration and index_line + 1 < len(file_lines):
                    index_line += 1
                    next_line = file_lines[index_line].strip()
                    if is_comment_only_line(next_line):
                        continue
                    if next_line == '};':
                        index_line -= 1
                        break
                    declaration += ' ' + next_line

                words = declaration.split()
                if len(words)>1 and words[0] in allowed_types:
                    variable_name = words[1].strip(';')
                    if len(words) > 2 and words[2].startswith('['):
                        variable_name += words[2].strip(';')
                    if '*' in variable_name:
                        # It is a pointer, like this: double * a_pointer.
                        index_line += 1
                        continue
                    structs.append('        ' + words[0] + ' ' + variable_name)
                elif len(words)>2 and words[0] == 'enum' and words[1] in enum_names:
                    variable_name = words[2].strip(';')
                    if len(words) > 3 and words[3].startswith('['):
                        variable_name += words[3].strip(';')
                    structs.append('        ' + words[1] + ' ' + variable_name)
                index_line += 1




# In[6]:

classes = []

class_names = ['FileContent', 'NonColdDarkMatter', 'InputModule', 'BackgroundModule',
                'ThermodynamicsModule', 'PerturbationsModule',
                'PrimordialModule', 'NonlinearModule', 'TransferModule', 'SpectraModule', 'LensingModule',
                'ClassConstants', 'NcdmSettings',]
allowed_types = ['double', 'int', 'short', 'char', 'bool', 'void', 'ErrorMsg', 'std::string',
                 'std::map<std::string, std::vector<double>>',
                 'std::map<std::string, int>', 'std::shared_ptr<NonColdDarkMatter>',
                 'std::vector<std::vector<double>>', 'std::vector<std::vector<short>>',
                 'std::vector<std::vector<int>>',
                 'std::vector<double>', 'std::vector<int>', 'std::vector<short>',
                 'std::vector<std::string>'] + struct_names

keywords_to_be_ignored = ['static', 'constexpr', 'const']

for file in h_files:
    file_lines = read_header_lines(file)
    class_name = ''
    index_line = 0
    while index_line < len(file_lines):
        line = file_lines[index_line]
        if not class_name:
            # Check for struct
            if 'class' not in line and 'struct' not in line:
                index_line += 1
                continue
            for s in class_names:
                if (('struct ' + s + ' ' in line or 'struct ' + s + '\n' in line) or
                    ('class ' + s + ' ' in line or 'class ' + s + '\n' in line)):
                    class_name = s
                    classes.append('cdef extern from "' + os.path.basename(file) + '":')
                    classes.append('    cdef cppclass ' + s + ':')
                    break
            index_line += 1
            continue

        if 'private:' in line or '};' in line:
            # Class just ended or we have reached the private section of the class
            class_name = ''
            classes.append('')
            index_line += 1
            continue

        if line.strip() in ['public:', 'protected:']:
            index_line += 1
            continue

        if line.strip() == '}':
            index_line += 1
            continue

        if is_comment_only_line(line):
            index_line += 1
            continue

        declaration = line.strip()
        while declaration and ';' not in declaration and '{' not in declaration and index_line + 1 < len(file_lines):
            index_line += 1
            next_line = file_lines[index_line].strip()
            if is_comment_only_line(next_line):
                continue
            if next_line == '}':
                continue
            if next_line in ['public:', 'private:', 'protected:', '};']:
                index_line -= 1
                break
            declaration += ' ' + next_line

        line = declaration
        # Cython does not like keywords
        removed_keyword = True
        while removed_keyword:
            removed_keyword = False
            for keyword in keywords_to_be_ignored:
                if line.startswith(keyword):
                    line = line[len(keyword):].strip()
                    removed_keyword = True
        typename = ''
        for sometype in allowed_types:
            if line.startswith(sometype):
                typename = sometype
                break
        if not typename:
            index_line += 1
            continue

        variable_name_begin = -1
        for index in range(len(typename), len(line)):
            if line[index] == '*' or line[index] == '&':
                typename += line[index]
            elif line[index] != ' ':
                variable_name_begin = index
                break
        if variable_name_begin == -1:
            index_line += 1
            continue

        variable_name_end = -1
        parantheses_count = 0
        bracket_count = 0
        is_function = False
        for index in range(variable_name_begin, len(line)):
            if ((line[index] == ' ' and parantheses_count == 0) or
                (line[index] == ';' and parantheses_count == 0)):
                if line[index] == ' ' and index + 1 < len(line) and line[index + 1] == '[':
                    continue
                if bracket_count > 0:
                    continue
                variable_name_end = index
                break
            if line[index] == '(':
                parantheses_count += 1
                is_function = True
            elif line[index] == ')':
                parantheses_count -= 1
            elif line[index] == '[':
                bracket_count += 1
            elif line[index] == ']':
                bracket_count -= 1

        if variable_name_end == -1:
            index_line += 1
            continue

        variable_name = line[variable_name_begin:variable_name_end]
        variable_name = variable_name.replace('enum ','')
        variable_name = variable_name.replace(' [', '[')
        if is_function:
            variable_name += ' except +'

        out_line = '        ' + typename + ' ' + variable_name
        out_line = out_line.replace('std::','').replace('<', '[').replace('>',']')
        classes.append(out_line)
        index_line += 1


# In[7]:


modules = [m for m in class_names if 'Module' in m]
cosmology_class = []
cosmology_class.append('cdef extern from "cosmology.h":')
for m in modules:
    cosmology_class.append('    ctypedef shared_ptr[const ' + m + '] ' + m + 'Ptr')
cosmology_class.append('')
cosmology_class.append('    cdef cppclass Cosmology:')
cosmology_class.append('        Cosmology(FileContent& fc) except +')
for m in modules:
    cosmology_class.append('        ' + m + 'Ptr& Get' + m + '()')


# In[8]:

with open(outdir / 'cclassy.pxd', 'w', encoding='utf-8') as fid:
    fid.write(preample)
    for lines in [enums, structs, classes]:
        # Replace defined variables:
        for index_line in range(len(lines)):
            for key, val in definitions_dict.items():
                lines[index_line] = lines[index_line].replace(key, val)
            # Cython spells C++ types without the std:: prefix and uses [] for
            # template brackets (e.g. std::string -> string). The classes section
            # already does this inline; apply it to all sections (incl. structs).
            lines[index_line] = lines[index_line].replace('std::', '').replace('<', '[').replace('>', ']')
        fid.write("\n".join(lines))
        fid.write("\n\n")
