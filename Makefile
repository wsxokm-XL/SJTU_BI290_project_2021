CC=gcc
CFLAGS=-Wall -std=gnu99

TARGET= 1.3
SRCS = 1.3.c

OBJS = $(SRCS:.c=.o)

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CC) -o $@ $^

.PHONY: clean
clean:
	rm -rf $(TARGET) $(OBJS)

%.o: %.c
	$(CC) $(CFLAGS) -o $@ -c $<

