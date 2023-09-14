
rate_io_bin = executable('test_rate_io','rate_io.cpp',
    include_directories: [common_include],
    dependencies: [json_dep],
    link_with: [common_lib]
)

# message(dependency('nlohmann_json').name())

test('rate_io', rate_io_bin)