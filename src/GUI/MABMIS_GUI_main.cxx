//#if defined(_MSC_VER)
//#pragma warning ( disable : 4786 )
//#pragma warning ( disable : 4996 ) //warning C4996: 'strcpy': This function or variable may be unsafe.
//// #pragma warning ( disable : 4244 )
//#endif

// c, c++
#if defined(_WIN32) && !defined(__CYGWIN__)
#include <windows.h>
#else
#include <unistd.h>
#endif

#include <iostream>
#include <string>
using namespace std;


// kwwidgets 
#include "vtkKWApplication.h"
#include "vtkKWWindowBase.h"
// vtk
#include <vtksys/SystemTools.hxx>
#include <vtksys/CommandLineArguments.hxx>
// my own
#include "vtkKWMABMISGUI.h"

bool isWin32 = false;
char curFolder[MAX_FILE_NAME_LENGTH];
char programFolder[MAX_FILE_NAME_LENGTH];

//extern "C" int KW_MABMIS_GUI_Lib_Init(Tcl_Interp *interp);

extern "C" int Kw_mabmis_gui_lib_Init(Tcl_Interp *interp);

int my_main(int argc, char *argv[])
{
	// Initialize Tcl
	Tcl_Interp *interp = vtkKWApplication::InitializeTcl(1, &argv[0], &cerr);
	if (!interp)
	{
		cerr << "Error: InitializeTcl failed" << endl ;
		return 1;
	}

	Kw_mabmis_gui_lib_Init(interp);

	// Create the application
	vtkKWApplication *app = vtkKWApplication::New();
	app->SetName("MABMIS_GUI");
	//app->RestoreApplicationSettingsFromRegistry();

	// Set a help link. Can be a remote link (URL), or a local file
	// vtksys::SystemTools::GetFilenamePath(__FILE__) + "/help.html";
	app->SetHelpDialogStartingPage("http://bric.unc.edu/ideagroup");

	// Add a window
	// Set 'SupportHelp' to automatically add a menu entry for the help link
	vtkKWWindowBase* window = vtkKWWindowBase::New();									//main window
	window->SupportHelpOn();
	app->AddWindow(window);
	//app->DisplayExitDialog();
	app->PromptBeforeExitOff();
	window->Create();
	window->SetSize(800,300);
	//window->SetDisplayPositionToScreenCenter();
	window->SetDisplayPositionToMasterWindowCenterFirst();
	// window->DisplayCloseDialog();
	// window->SetPromptBeforeClose();
	// window->PromptBeforeCloseOn();

	vtkKWMABMISGUI* MABMISGUI = vtkKWMABMISGUI::New();
	MABMISGUI->SetApplication(app);
	MABMISGUI->SetWin(window);
	MABMISGUI->ShowHelpAboutDialog();
	//MABMISGUI->GlobalWarningDisplayOff();
#if defined(_WIN32) && !defined(__CYGWIN__)	// win32
	// to get current folder
	//HMODULE module = GetModuleHandle(0); 
	//CHAR buf[MAX_FILE_NAME_LENGTH]; 
	//GetModuleFileName(module, buf, sizeof buf); 
	//strcpy(programFolder, buf);
	//std::cerr << "programFolder: " << programFolder << std::endl;

	//TCHAR szDirectory[MAX_FILE_NAME_LENGTH] = "";
	//if(!::GetCurrentDirectory(sizeof(curFolder) - 1, curFolder))
	DWORD dwRet;
	dwRet = GetCurrentDirectory(sizeof(curFolder) - 1, curFolder);
	//std::cerr << "curFolder: " << curFolder << std::endl;
	//strcpy(curFolder,szDirectory);


	// to get programFolder
	//char* curDir = new char[MAX_FILE_NAME_LENGTH];
	//curDir[0] = '\0';
	const char* path = app->GetInstallationDirectory();
	//std::cerr << "path: " << path << std::endl;
	string path1 = vtksys::SystemTools::GetRealPath(path);
	//std::cerr << "path1: " << path1 <<std::endl;
	//MABMISGUI->SetCurrentDirectory(path1.c_str());
	strcpy(programFolder, path1.c_str());
	//std::cerr << "programFolder: " << programFolder << std::endl;


#else //unix64
	// to get current folder
	getcwd(curFolder, MAX_FILE_NAME_LENGTH);
	//std::cerr << "curFolder: " << curFolder << std::endl;

	// to get programFolder
	const char* path = app->GetInstallationDirectory();
	//std::cerr << "path: " << path << std::endl;
	string path1 = vtksys::SystemTools::GetRealPath(path);
	//std::cerr << "path1: " << path1 <<std::endl;
	//MABMISGUI->SetCurrentDirectory(path1.c_str());
	strcpy(programFolder, path1.c_str());
	//std::cerr << "programFolder: " << programFolder << std::endl;

#endif
	MABMISGUI->SetCurrentDirectory(curFolder);

	MABMISGUI->Run();

	//start application
	window->Display();
	app->Start();
	window->Close();

	MABMISGUI->Delete();
	
	window->Delete();
	app->Delete();

	Tcl_DeleteInterp(interp);
	Tcl_Finalize();
	
	//delete curDir;

	return 0;

}

#if defined(_WIN32) && !defined(__CYGWIN__)
#include <windows.h>
int __stdcall WinMain(HINSTANCE, HINSTANCE, LPSTR lpCmdLine, int)
{
	int argc;
	char **argv;
	vtksys::SystemTools::ConvertWindowsCommandLineToUnixArguments(lpCmdLine, &argc, &argv);
	isWin32 = true;
	int ret = my_main(argc, argv);
	for (int i = 0; i < argc; i++) { delete [] argv[i]; }
	delete [] argv;
	return ret;
}
#else
int main(int argc, char *argv[])
{
	isWin32 = false;
	return my_main(argc, argv);
}
#endif
