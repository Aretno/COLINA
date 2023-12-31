# -*-Makefile-*-

# Needed to make `"$@"` usable in recipes
set positional-arguments := true

default:
  just run --beam-on 10

test *REST:
  meson setup build/test test
  meson compile -C build/test
  meson install -C build/test
  sh install/test/run-each-test-in-separate-process.sh "$@"

catch2-demo *REST:
  echo "$@"
  meson setup build/test test
  meson compile -C build/test
  meson install -C build/test
  install/test/bin/catch2-demo-test "$@"

build:
  meson setup build/app src
  meson compile -C build/app

install: build
  meson install -C build/app

run *ARGS: install
  #!/usr/bin/env sh
  sh execute-with-nixgl-if-needed.sh ./install/app/bin/COLINA "$@"
  exit $?

clean:
  rm build install -rf
