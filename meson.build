project(
  'cMag',
  'c',
  version: '0.1.0',
  default_options: [
    'c_std=c99',
    'warning_level=2',
  ],
)

# Match your Makefile's -D_XOPEN_SOURCE=600
add_project_arguments('-D_XOPEN_SOURCE=600', language: 'c')

inc = include_directories('includes')

lib_sources = files(
  'src/maggrid.c',
  'src/magfield.c',
  'src/magfieldutil.c',
  'src/mapcolor.c',
  'src/magfielddraw.c',
  'src/magfieldio.c',
  'src/svg.c',
  'src/testdata.c',
)

lib_cmag = static_library(
  'cMag',                    # -> libcMag.a
  lib_sources,
  include_directories: inc,
  install: true,
)

cc = meson.get_compiler('c')
m_dep = cc.find_library('m', required: false)  # -lm (Linux); no-op on macOS

executable(
  'cMagTest',
  'src/main.c',
  include_directories: inc,
  link_with: lib_cmag,
  dependencies: m_dep,
  install: true,
)

# Optional: install headers (recommended)
install_headers(
  'includes/magfield.h',
  'includes/magfielddraw.h',
  'includes/magfieldio.h',
  'includes/magfieldutil.h',
  'includes/maggrid.h',
  'includes/mapcolor.h',
  'includes/munittest.h',
  'includes/svg.h',
  'includes/testdata.h',
)

