
#include "MyClassCent.C"
#include "READCent.C"

void start(){
MyClass t("/home/dim2/FLOW5/chain/chainBIG.root");
t.Book();t.Loop();
t.SaveData("~/FLOW5/OUT/urqmdHome.root");
}
void macroLoopH(){start();
read("/home/dim2/FLOW5/OUT/urqmdHome.root"); 
}
