
#include "MyClass.C"
#include "READ.C"

void start(){
MyClass t("/home/dim2/FLOW5/All.root");
t.Book();t.Loop();
t.SaveData("~/FLOW5/OUT/urqmd.root");
}
void macroLoop(){//start();
read("~/FLOW5/OUT/urqmd.root"); }
