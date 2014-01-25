all:
	g++ -O3 -I /usr/local/include/ -L /usr/local/lib -o sc2d.x -framework OpenGL -framework Cocoa -framework IOKit sc2d.cpp -lglfw

clean:
	rm *.x
