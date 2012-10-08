#ifndef __vtkKWMABMISGUI_h
#define __vtkKWMABMISGUI_h

#include "vtkKWObject.h"

class vtkKWWindow;
class vtkKWWindowBase;


class vtkKWPushButtonWithLabel;
class vtkKWCheckButtonWithLabel;
class vtkKWEntryWithLabel;
class vtkKWSpinBox;
class vtkKWSpinBoxWithLabel;
class vtkKWComboBox;
class vtkKWComboBoxWithLabel;
class vtkKWComboBoxSet;
class vtkKWFrameWithLabel;
class vtkKWFileBrowserWidget;
class vtkKWDirectoryExplorer;
class vtkKWFileBrowserDialog;
class vtkKWFileListTable;
class vtkKWLabel;
class vtkKWRadioButton;
class vtkKWRadioButtonSet;


#include <iostream>
#include <string>
using namespace std;
#define MAX_FILE_NAME_LENGTH 1024
extern bool isWin32;
extern char curFolder[];
extern char programFolder[];

class vtkKWMABMISGUI: public vtkKWObject
{

public:

	vtkKWMABMISGUI();
	~vtkKWMABMISGUI();

	static vtkKWMABMISGUI* New();
	vtkTypeRevisionMacro(vtkKWMABMISGUI,vtkKWObject);

	virtual int Run();
	virtual void ShowHelpAboutDialog();
	void SetWin(vtkKWWindowBase *win)
	{window = win;}
	void SetCurrentDirectory(char* dir)
	{
		strcpy(m_cCurrentDirectory, dir);
	}
	//variable
protected:

	vtkKWWindowBase *window;

	vtkKWFrameWithLabel*				m_flPara;
	vtkKWFrameWithLabel*				m_flData;
	vtkKWFrameWithLabel*				m_flInfo;

	/////////////////////////////////////
	vtkKWPushButtonWithLabel*			m_pbStartMABMIS;
	//
	//vtkKWSpinBoxWithLabel*				m_sbRunGroupMean;
	////
	//vtkKWCheckButtonWithLabel*			m_cbHistogramMatch;
	////
	//vtkKWCheckButtonWithLabel*			m_cbAffineRegistration;
	////
	//vtkKWSpinBoxWithLabel*				m_sbLevelMax;
	////
	////vtkKWFrameWithLabel*				m_flSmoothKernel;
	//vtkKWLabel*							m_lSmoothKernel;
	////vtkKWComboBox*						m_cbSigmaInit;
	////vtkKWComboBox*						m_cbSigmaMid;
	////vtkKWComboBox*						m_cbSigmaLast;
	////vtkKWLabel*							m_lSigmaInit;
	////vtkKWLabel*							m_lSigmaMid;
	////vtkKWLabel*							m_lSigmaLast;
	//vtkKWComboBoxWithLabel*					m_cblSigmaInit;
	//vtkKWComboBoxWithLabel*					m_cblSigmaMid;
	//vtkKWComboBoxWithLabel*					m_cblSigmaLast;

	vtkKWEntryWithLabel*				m_elAtlasesSize;
	//
	//vtkKWSpinBoxWithLabel*				m_sbNeighborhoodSize;
	////
	//vtkKWRadioButtonSet*				m_rbsOutputType;
	//
	//vtkKWLabel*							m_lChooseDirectory;
	vtkKWLabel*							m_lWarningInfo;
	vtkKWFileBrowserWidget*				m_fbwChooseDirectory;
	vtkKWLabel*							m_lBlankInfo;

	bool								m_bFilesValid;
	char*								m_cCurrentDirectory;
	
public:

	virtual void MABMISRunCallback();

protected:

	void InitializePara();
	void CheckFileValidity();
};


#endif