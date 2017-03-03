CFLAGS=-Wall -02 -g -Wextra -Isrc -rdynamic - DNDEBUG $(OPTFLAGS)

SOURCES=$(wildcard src/**/*.c src/*.c)
OBJECTS=$(patsubst %.c,%.o,$(TEST_SRC))

TARGET=build/lib_alg_lin.a
SO_TARGET=$(patsubst %.a,%.so,$(TARGET))

## The target Build
all: $(TARGET) $(SO_TARGET) 

dev: CFLAGS=-g -Wall -Isrc -Wextra $(OPTFLAGS)
dev: all

$(TARGET): CFLAGS += -fPIC
$(TARGET): build $(OBJECTS)
	ar rcs $@ $(OBJECTS)
	ranlib $@

$(SO_TARGET): $(TARGET) $(OBJECTS)
	$(CC) -shared -o $@ $(OBEJCTS)

build:
	@mkdir -p build
	@mkdir -p bin

## The Cleaner
clean:
	rm -rf build $(OBJECTS) $(TESTS)
	find . -name "*.gc*" -exec rm {} \;
	rm -rf 'find . -name "*.dSYM" -print'

## The Install
install: all
	install -d $(DESTDIR)/$(PREFIX)/lib
	install $(TARGET) $(DESTDIR)/$(PREFIX)/lib

## The checker
 ##
