beng207_2026_sources = [
  'main.c',
  'beng207_2026-application.c',
  'beng207_2026-window.c',
]

beng207_2026_deps = [
  dependency('gtk4'),
  dependency('libadwaita-1', version: '>= 1.4'),
]

beng207_2026_sources += gnome.compile_resources('beng207_2026-resources',
  'beng207_2026.gresource.xml',
  c_name: 'beng207_2026'
)

cc = meson.get_compiler('c')
m_dep = cc.find_library('m', required: true)
beng207_2026_deps += [m_dep]

executable('beng207_2026', beng207_2026_sources,
  dependencies: beng207_2026_deps,
  install: true,
)
