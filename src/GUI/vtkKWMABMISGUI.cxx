// update 20100714: GUI update with correct expansion and arrangement
// update 20100729: GUI neighbor size change from 3 to 5 (previously 1 to 5, which is not consistent with command version.)
#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#pragma warning ( disable : 4996 ) //warning C4996: 'strcpy': This function or variable may be unsafe.
// #pragma warning ( disable : 4244 )
#endif

#include "vtkKWMABMISGUI.h"

//head file of KWWidgets
#include "vtkObjectFactory.h"
#include "vtkKWApplication.h"
#include "vtkKWWindowBase.h"

#include "vtkKWPushButton.h"
#include "vtkKWPushButtonWithLabel.h"
#include "vtkKWCheckButton.h"
#include "vtkKWCheckButtonWithLabel.h"
#include "vtkKWEntry.h"
#include "vtkKWEntryWithLabel.h"
#include "vtkKWSpinBox.h"
#include "vtkKWSpinBoxWithLabel.h"
#include "vtkKWComboBox.h"
#include "vtkKWComboBoxWithLabel.h"
#include "vtkKWComboBoxSet.h"
#include "vtkKWFrame.h"
#include "vtkKWFrameWithLabel.h"
#include "vtkKWFileBrowserWidget.h"
#include "vtkKWDirectoryExplorer.h"
#include "vtkKWFileBrowserDialog.h"
#include "vtkKWFileListTable.h"
#include "vtkKWLabel.h"
#include "vtkKWRadioButton.h"
#include "vtkKWRadioButtonSet.h"

#include "vtkKWMenu.h"
#include "vtkKWMessageDialog.h"

#include "vtkToolkits.h"
#include <vtksys/SystemTools.hxx>


//c++ lib
#include <iostream>
#include <sstream>
#include <fstream>
#include <string.h>


vtkStandardNewMacro( vtkKWMABMISGUI );
vtkCxxRevisionMacro( vtkKWMABMISGUI, "$Revision: 1.0$");
vtkKWMABMISGUI::vtkKWMABMISGUI()
{
	
	m_flPara	= NULL;
	m_flData	= NULL;
	m_flInfo	= NULL;

	////////////////////////////
	m_pbStartMABMIS			= NULL;
	//
	//m_sbRunGroupMean		= NULL;
	////
	//m_cbHistogramMatch		= NULL;
	////
	//m_cbAffineRegistration	= NULL;
	////
	//m_sbLevelMax			= NULL;
	////
	////m_flSmoothKernel		= NULL;
	//m_lSmoothKernel			= NULL;

	//m_cblSigmaInit			= NULL;
	//m_cblSigmaMid			= NULL;
	//m_cblSigmaLast			= NULL;
	////m_lSigmaInit			= NULL;
	////m_lSigmaMid				= NULL;
	////m_lSigmaLast			= NULL;
	m_elAtlasesSize			= NULL;
	//
	//m_sbNeighborhoodSize	= NULL;
	////
	//m_rbsOutputType			= NULL;
	//
	m_fbwChooseDirectory	= NULL;
	//m_lChooseDirectory		= NULL;
	m_lWarningInfo			= NULL;

	m_lBlankInfo			= NULL;

	m_cCurrentDirectory		= NULL;
	InitializePara();
}
vtkKWMABMISGUI::~vtkKWMABMISGUI()
{
	if(this->m_flPara)			this->m_flPara->Delete();
	if(this->m_flData)			this->m_flData->Delete();
	if(this->m_flInfo)			this->m_flInfo->Delete();


	///////////////////////////////////////
	if(this->m_pbStartMABMIS)			this->m_pbStartMABMIS->Delete();
	//
	//if(this->m_sbRunGroupMean)			this->m_sbRunGroupMean->Delete();
	////
	//if(this->m_cbHistogramMatch)		this->m_cbHistogramMatch->Delete();
	////
	//if(this->m_cbAffineRegistration)	this->m_cbAffineRegistration->Delete();
	////
	//if(this->m_sbLevelMax)				this->m_sbLevelMax->Delete();
	////
	////if(this->m_flSmoothKernel)		this->m_flSmoothKernel->Delete();
	//if(this->m_lSmoothKernel)			this->m_lSmoothKernel->Delete();
	//if(this->m_cblSigmaInit)			this->m_cblSigmaInit->Delete();
	//if(this->m_cblSigmaMid)				this->m_cblSigmaMid->Delete();
	//if(this->m_cblSigmaLast)			this->m_cblSigmaLast->Delete();
	////if(this->m_lSigmaInit)			this->m_lSigmaInit->Delete();
	////if(this->m_lSigmaMid)				this->m_lSigmaMid->Delete();
	////if(this->m_lSigmaLast)			this->m_lSigmaLast->Delete();
	if(this->m_elAtlasesSize)			this->m_elAtlasesSize->Delete();
	//
	//if(this->m_sbNeighborhoodSize)		this->m_sbNeighborhoodSize->Delete();
	//
	//if(this->m_rbsOutputType)			this->m_rbsOutputType->Delete();
	//
	if(this->m_fbwChooseDirectory)		this->m_fbwChooseDirectory->Delete();
	//if(this->m_lChooseDirectory)		this->m_lChooseDirectory->Delete();
	if(this->m_lWarningInfo)			this->m_lWarningInfo->Delete();
	if(this->m_lBlankInfo)				this->m_lBlankInfo->Delete();

	if(this->m_cCurrentDirectory)			delete this->m_cCurrentDirectory;
}

int vtkKWMABMISGUI::Run()
{
	double readonlycolor[3] = {1.0,1.0,1.0};
	int window_width = this->window->GetWidth();
	int window_height = this->window->GetHeight();

	vtkKWApplication *app = this->GetApplication();
	vtkKWFrame* frame_main = this->window->GetViewFrame();

	vtkKWMenu * help_menu = this->window->GetHelpMenu();
	help_menu->SetItemCommand(help_menu->GetNumberOfItems()-1, this, "ShowHelpAboutDialog");
	
	/////////////////////////////////////////////////////////////////
	// set Info frame
	if (!m_flInfo)
		m_flInfo = vtkKWFrameWithLabel::New();
	m_flInfo->SetParent(frame_main);
	m_flInfo->Create();
	m_flInfo->SetLabelText("System information");
	m_flInfo->GetLabel()->SetFont("courier 10 normal");
	m_flInfo->AllowFrameToCollapseOff();
	m_flInfo->SetDefaultLabelFontWeightToNormal();
	m_flInfo->SetBalloonHelpString("Show system information about MABMIS.");
	m_flInfo->SetHeight(40);
	m_flInfo->SetWidth(800);
	//app->Script("place %s -width 200 -height 120 -x 20 -y 70", m_flPara->GetWidgetName());

	// set Para frame
	if (!m_flPara)
		m_flPara = vtkKWFrameWithLabel::New();
	m_flPara->SetParent(frame_main);
	m_flPara->Create();
	m_flPara->SetLabelText("Choose parameters");
	m_flPara->GetLabel()->SetFont("courier 10 normal");
	m_flPara->AllowFrameToCollapseOff();
	m_flPara->SetDefaultLabelFontWeightToNormal();
	m_flPara->SetBalloonHelpString("Choose all parameters in MABMIS.");
	m_flPara->SetHeight(200);
	m_flPara->SetWidth(250);
	//m_flPara->GetFrame()->SetHeight(400);
	//app->Script("place %s -width 200 -height 120 -x 20 -y 70", m_flPara->GetWidgetName());
	//origin
	//app->Script("pack %s -side left -anchor nw -expand y -padx 2 -pady 2", m_flPara->GetWidgetName());

	// set Data frame
	if (!m_flData)
		m_flData = vtkKWFrameWithLabel::New();
	m_flData->SetParent(frame_main);
	m_flData->Create();
	m_flData->SetLabelText("Choose image ID file (.txt)");
	m_flData->GetLabel()->SetFont("courier 10 normal");
	m_flData->AllowFrameToCollapseOff();
	m_flData->SetDefaultLabelFontWeightToNormal();
	m_flData->SetBalloonHelpString("Choose image ID file (.txt)");
	m_flData->SetHeight(200);
	m_flData->SetWidth(550);
	//app->Script("place %s -width 200 -height 120 -x 20 -y 70", m_flPara->GetWidgetName());
	//original
	//app->Script("pack %s -side left -anchor nw -expand y -padx 2 -pady 2", m_flData->GetWidgetName());

	//app->Script("pack %s -side top  -anchor nw -expand 0 -fill x -padx 2 -pady 2", m_flInfo->GetWidgetName());
	//app->Script("pack %s -side left -anchor c  -expand 1 -fill y -padx 2 -pady 2", m_flPara->GetWidgetName());
	//app->Script("pack %s -side left -anchor c  -expand 1 -fill both -padx 2 -pady 2", m_flData->GetWidgetName());

	app->Script("pack %s -side top  -anchor nw -expand 0 -fill x -padx 2 -pady 2", m_flInfo->GetWidgetName());
	app->Script("pack %s -side left -anchor c  -expand 0 -fill y -padx 2 -pady 2", m_flPara->GetWidgetName());
	app->Script("pack %s -side left -anchor c  -expand 1 -fill both -padx 2 -pady 2", m_flData->GetWidgetName());

	///////////////////////////////////////////////

	if (!m_lWarningInfo)
		m_lWarningInfo = vtkKWLabel::New();
	m_lWarningInfo->SetParent(m_flInfo->GetFrame());
	m_lWarningInfo->Create();
	m_lWarningInfo->SetText("Select image ID lists file (.txt)!");
	m_lWarningInfo->SetForegroundColor(0,0,0);
	m_lWarningInfo->SetFont("system 12 normal");
	m_lWarningInfo->SetWidth(600);
	app->Script("pack %s -side top -anchor c -expand y -padx 2 -pady 2 -in %s", 
		m_lWarningInfo->GetWidgetName(), m_flInfo->GetFrame()->GetWidgetName());

	/////////////////////////////////////////////////////////////////
	// blank
	if (!m_lBlankInfo)
		m_lBlankInfo = vtkKWLabel::New();
	m_lBlankInfo->SetParent(m_flPara->GetFrame());
	m_lBlankInfo->Create();
	m_lBlankInfo->SetText("Set the number of atlases:");
	m_lBlankInfo->SetFont("courier 12 normal");
	//m_lBlankInfo->SetWidth(10);
	m_lBlankInfo->SetJustificationToLeft();
	m_lBlankInfo->SetBalloonHelpString("Set the number of atlases:");
	app->Script("pack %s -side top -anchor nw -expand 1 -padx 5 -pady 0 -in %s", 
		m_lBlankInfo->GetWidgetName(), m_flPara->GetFrame()->GetWidgetName());

	// set widget for atlases size
	if (!m_elAtlasesSize)
		m_elAtlasesSize = vtkKWEntryWithLabel::New();
	m_elAtlasesSize->SetParent(m_flPara->GetFrame());
	m_elAtlasesSize->Create();
	m_elAtlasesSize->GetWidget()->SetRestrictValueToInteger();
	m_elAtlasesSize->GetWidget()->SetValueAsInt(10);
	m_elAtlasesSize->GetWidget()->SetWidth(5);
	m_elAtlasesSize->GetWidget()->SelectAllOnFocusInOn();
	m_elAtlasesSize->SetBalloonHelpString("Choose the number of atlases in all images. ");
	m_elAtlasesSize->SetLabelText("Number of atlases: ");
	m_elAtlasesSize->GetLabel()->SetFont("Courier 12 normal");
	//app->Script("place %s -width 180 -height 20 -x 20 -y 20", m_sbNeighborhoodSize->GetWidgetName());
	app->Script("pack %s -side top -anchor nw -expand 1 -padx 5 -pady 5 -in %s", 
		m_elAtlasesSize->GetWidgetName(), m_flPara->GetFrame()->GetWidgetName());


/*	// set widget for neighborhood size
	if(!m_sbNeighborhoodSize)
		m_sbNeighborhoodSize = vtkKWSpinBoxWithLabel::New();
	m_sbNeighborhoodSize->SetParent(m_flPara->GetFrame());
	m_sbNeighborhoodSize->Create();
	m_sbNeighborhoodSize->GetWidget()->SetRestrictValueToInteger();
	m_sbNeighborhoodSize->GetWidget()->SetRange(3, 5);
	m_sbNeighborhoodSize->GetWidget()->SetValue(3);
	m_sbNeighborhoodSize->GetWidget()->SetIncrement(1);
	m_sbNeighborhoodSize->GetWidget()->SetRestrictValueToInteger();
	m_sbNeighborhoodSize->GetWidget()->SetWidth(6);
	m_sbNeighborhoodSize->GetWidget()->SetStateToReadOnly();
	m_sbNeighborhoodSize->GetWidget()->SetReadOnlyBackgroundColor(readonlycolor);//
	m_sbNeighborhoodSize->SetBalloonHelpString("Choose the number of nearest neighbors. ");
	m_sbNeighborhoodSize->SetLabelText("Neighborhood size: ");
	m_sbNeighborhoodSize->GetLabel()->SetFont("Courier 12 normal");
	//app->Script("place %s -width 180 -height 20 -x 20 -y 20", m_sbNeighborhoodSize->GetWidgetName());
	app->Script("pack %s -side top -anchor nw -expand 1 -padx 5 -pady 5 -in %s", 
		m_sbNeighborhoodSize->GetWidgetName(), m_flPara->GetFrame()->GetWidgetName());


	// set Smooth Kernel frame
	if (!m_lSmoothKernel)
		m_lSmoothKernel = vtkKWLabel::New();
	m_lSmoothKernel->SetParent(m_flPara->GetFrame());
	m_lSmoothKernel->Create();
	m_lSmoothKernel->SetText("Size of smooth kernels");
	m_lSmoothKernel->SetFont("courier 12 italic");
	m_lSmoothKernel->SetJustificationToLeft();
	m_lSmoothKernel->SetBalloonHelpString("Choose different smooth kernel in Demons. ");
	//app->Script("place %s -width 200 -height 120 -x 20 -y 70", m_lSmoothKernel->GetWidgetName());
	app->Script("pack %s -side top -anchor nw -expand 1 -padx 5 -pady 2 -in %s", 
		m_lSmoothKernel->GetWidgetName(), m_flPara->GetFrame()->GetWidgetName());

	const char* sigma_value[] = {"1.0", "1.5", "2.0", "2.5", "3.0"};
	if (!m_cblSigmaInit)
		m_cblSigmaInit = vtkKWComboBoxWithLabel::New();
	m_cblSigmaInit->SetParent(m_flPara->GetFrame());
	m_cblSigmaInit->Create();
	for (int i = 0; i < sizeof(sigma_value) / sizeof(sigma_value[0]); i++)
	{
		m_cblSigmaInit->GetWidget()->AddValue(sigma_value[i]);
	}
	m_cblSigmaInit->GetWidget()->ReadOnlyOn();
	m_cblSigmaInit->GetWidget()->SetValue("2.0");
	m_cblSigmaInit->GetWidget()->SetWidth(5);
	m_cblSigmaInit->SetLabelText("    Initial level: ");
	m_cblSigmaInit->GetLabel()->SetFont("courier 12 normal");
	//m_cblSigmaInit->SetWidth(100);
	//app->Script("place %s -width 80 -height 20 -x 120 -y 90", m_cbSigmaInit->GetWidgetName());
	app->Script("pack %s -side top -anchor nw -expand 1 -padx 2 -pady 2 -in %s", 
		m_cblSigmaInit->GetWidgetName(), m_flPara->GetFrame()->GetWidgetName());


	if (!m_cblSigmaMid)
		m_cblSigmaMid = vtkKWComboBoxWithLabel::New();
	m_cblSigmaMid->SetParent(m_flPara->GetFrame());
	m_cblSigmaMid->Create();
	for (int i = 0; i < sizeof(sigma_value) / sizeof(sigma_value[0]); i++)
	{
		m_cblSigmaMid->GetWidget()->AddValue(sigma_value[i]);
	}
	m_cblSigmaMid->GetWidget()->ReadOnlyOn();
	m_cblSigmaMid->GetWidget()->SetValue("2.0");
	m_cblSigmaMid->GetWidget()->SetWidth(5);
	m_cblSigmaMid->SetLabelText("    Middle levels: ");
	m_cblSigmaMid->GetLabel()->SetFont("courier 12 normal");
	//m_cblSigmaMid->SetWidth(100);
	//app->Script("place %s -width 80 -height 20 -x 120 -y 120", m_cblSigmaMid->GetWidgetName());
	app->Script("pack %s -side top -anchor nw -expand 1 -padx 2 -pady 2 -in %s", 
		m_cblSigmaMid->GetWidgetName(), m_flPara->GetFrame()->GetWidgetName());


	if (!m_cblSigmaLast)
		m_cblSigmaLast = vtkKWComboBoxWithLabel::New();
	m_cblSigmaLast->SetParent(m_flPara->GetFrame());
	m_cblSigmaLast->Create();
	for (int i = 0; i < sizeof(sigma_value) / sizeof(sigma_value[0]); i++)
	{
		m_cblSigmaLast->GetWidget()->AddValue(sigma_value[i]);
	}
	m_cblSigmaLast->GetWidget()->ReadOnlyOn();
	m_cblSigmaLast->GetWidget()->SetValue("2.0");
	m_cblSigmaLast->GetWidget()->SetWidth(5);
	m_cblSigmaLast->SetLabelText("    Last level:    ");
	m_cblSigmaLast->GetLabel()->SetFont("courier 12 normal");
	//m_cblSigmaLast->SetWidth(100);

	//app->Script("place %s -width 80 -height 20 -x 120 -y 150", m_cblSigmaLast->GetWidgetName());
	app->Script("pack %s -side top -anchor nw -expand 1 -padx 2 -pady 2 -in %s", 
		m_cblSigmaLast->GetWidgetName(), m_flPara->GetFrame()->GetWidgetName());

	// set maximum levels
	if(!m_sbLevelMax)
		m_sbLevelMax = vtkKWSpinBoxWithLabel::New();
	m_sbLevelMax->SetParent(m_flPara->GetFrame());
	m_sbLevelMax->Create();
	m_sbLevelMax->GetWidget()->SetRestrictValueToInteger();
	m_sbLevelMax->GetWidget()->SetRange(3, 5);
	m_sbLevelMax->GetWidget()->SetValue(3);
	m_sbLevelMax->GetWidget()->SetIncrement(1);
	m_sbLevelMax->GetWidget()->SetRestrictValueToInteger();
	m_sbLevelMax->GetWidget()->SetWidth(6);
	m_sbLevelMax->GetWidget()->SetStateToReadOnly();
	m_sbLevelMax->GetWidget()->SetReadOnlyBackgroundColor(readonlycolor);//
	m_sbLevelMax->SetBalloonHelpString("Choose the number of max levels. ");
	m_sbLevelMax->SetLabelText("Maximum levels:    ");
	m_sbLevelMax->GetLabel()->SetFont("courier 12 normal");
	//app->Script("place %s -width 180 -height 20 -x 20 -y 220", m_sbLevelMax->GetWidgetName());
	app->Script("pack %s -side top -anchor nw -expand 1 -padx 5 -pady 2 -in %s", 
		m_sbLevelMax->GetWidgetName(), m_flPara->GetFrame()->GetWidgetName());
	
	// set entry with label of additional registration to group mean
	if (!m_sbRunGroupMean)
		m_sbRunGroupMean = vtkKWSpinBoxWithLabel::New();
	m_sbRunGroupMean->SetParent(m_flPara->GetFrame());
	m_sbRunGroupMean->Create();
	m_sbRunGroupMean->GetWidget()->SetRestrictValueToInteger();
	m_sbRunGroupMean->GetWidget()->SetRange(0, 3);
	m_sbRunGroupMean->GetWidget()->SetValue(1);
	m_sbRunGroupMean->GetWidget()->SetIncrement(1);
	m_sbRunGroupMean->GetWidget()->SetRestrictValueToInteger();
	m_sbRunGroupMean->GetWidget()->SetWidth(2);
	m_sbRunGroupMean->GetWidget()->SetStateToReadOnly();
	m_sbRunGroupMean->GetWidget()->SetReadOnlyBackgroundColor(readonlycolor);//
	m_sbRunGroupMean->SetBalloonHelpString("Set the number of additional runs of registration to the mean image. ");
	m_sbRunGroupMean->SetLabelText("Registration to mean:");
	m_sbRunGroupMean->GetLabel()->SetFont("courier 12 normal");
	//app->Script("place %s -width 180 -height 20 -x 20 -y 220", m_sbLevelMax->GetWidgetName());
	app->Script("pack %s -side top -anchor nw -expand 1 -padx 5 -pady 5 -in %s", 
		m_sbRunGroupMean->GetWidgetName(), m_flPara->GetFrame()->GetWidgetName());

	// set check button for histogram matching
	if (!m_cbHistogramMatch)
		m_cbHistogramMatch = vtkKWCheckButtonWithLabel::New();
	m_cbHistogramMatch->SetParent(m_flPara->GetFrame());
	m_cbHistogramMatch->Create();
	m_cbHistogramMatch->SetLabelText("Histogram matching:  ");
	m_cbHistogramMatch->GetLabel()->SetFont("courier 12 normal");
	m_cbHistogramMatch->GetWidget()->SetSelectedState(1);
	m_cbHistogramMatch->SetBalloonHelpString("Check if do histogram matching.");
	//app->Script("place %s -width 180 -height 20 -x 20 -y 250", m_cbHistogramMatch->GetWidgetName());
	app->Script("pack %s -side top -anchor nw -expand 1 -padx 5 -pady 5 -in %s", 
		m_cbHistogramMatch->GetWidgetName(), m_flPara->GetFrame()->GetWidgetName());

	// set check button for affine registration
	if (!m_cbAffineRegistration)
		m_cbAffineRegistration = vtkKWCheckButtonWithLabel::New();
	m_cbAffineRegistration->SetParent(m_flPara->GetFrame());
	m_cbAffineRegistration->Create();
	m_cbAffineRegistration->SetLabelText("Affine registration: ");
	m_cbAffineRegistration->GetLabel()->SetFont("courier 12 normal");
	m_cbAffineRegistration->GetWidget()->SetSelectedState(0);
	m_cbAffineRegistration->SetBalloonHelpString("Check if do affine registration before non-rigid.");
	//app->Script("place %s -width 180 -height 20 -x 20 -y 250", m_cbHistogramMatch->GetWidgetName());
	app->Script("pack %s -side top -anchor nw -expand 1 -padx 5 -pady 5 -in %s", 
		m_cbAffineRegistration->GetWidgetName(), m_flPara->GetFrame()->GetWidgetName());
*/
/*
	// radio button set for output type
	if (!m_rbsOutputType)
		m_rbsOutputType = vtkKWRadioButtonSet::New();
	m_rbsOutputType->SetParent(m_flPara->GetFrame());
	m_rbsOutputType->Create();
	m_rbsOutputType->SetBorderWidth(2);
	m_rbsOutputType->SetReliefToGroove();
	m_rbsOutputType->SetBalloonHelpString("Choose the output data type.");
	m_rbsOutputType->SetMaximumNumberOfWidgetsInPackingDirection(2);

	char outputTypeFloat[]="float";
	char outputTypeInt[]="int";
	char outputTypeShort[]="short";
	char outputTypeUChar[]="unsigned char";

	vtkKWRadioButton *radiob0 = m_rbsOutputType->AddWidget(0);
	radiob0->SetText(outputTypeFloat);
	radiob0->SetFont("courier 12 normal");
	radiob0->SetBalloonHelpString("Choose the output data type.");
	vtkKWRadioButton *radiob1 = m_rbsOutputType->AddWidget(1);
	radiob1->SetText(outputTypeInt);
	radiob1->SetFont("courier 12 normal");
	radiob1->SetBalloonHelpString("Choose the output data type.");
	vtkKWRadioButton *radiob2 = m_rbsOutputType->AddWidget(2);
	radiob2->SetText(outputTypeShort);
	radiob2->SetFont("courier 12 normal");
	radiob2->SetBalloonHelpString("Choose the output data type.");
	vtkKWRadioButton *radiob3 = m_rbsOutputType->AddWidget(3);
	radiob3->SetText(outputTypeUChar);
	radiob3->SetFont("courier 12 normal");
	radiob3->SetBalloonHelpString("Choose the output data type.");

	m_rbsOutputType->GetWidget(3)->SetSelectedState(1);
	app->Script("pack %s -side top -anchor nw -expand 1 -padx 5 -pady 2 -in %s", 
		m_rbsOutputType->GetWidgetName(), m_flPara->GetFrame()->GetWidgetName());

*/
	// set widget for start button
	if(!m_pbStartMABMIS)
		m_pbStartMABMIS = vtkKWPushButtonWithLabel::New();
	m_pbStartMABMIS->SetParent(m_flPara->GetFrame());
	//m_pbStartMABMIS->SetParent(frame_main);
	m_pbStartMABMIS->Create();
	m_pbStartMABMIS->GetWidget()->SetText("Run");
	m_pbStartMABMIS->SetLabelText("Start MABMIS:   ");
	m_pbStartMABMIS->GetLabel()->SetFont("courier 12 normal");
	m_pbStartMABMIS->SetBalloonHelpString("Choose the image list file and press this button to run MABMIS. ");
	m_pbStartMABMIS->GetWidget()->SetCommand(this, "MABMISRunCallback");
	m_pbStartMABMIS->GetWidget()->SetWidth(10);
	app->Script("pack %s -side top -anchor nw -expand 1 -padx 5 -pady 5 -in %s", 
		m_pbStartMABMIS->GetWidgetName(), m_flPara->GetFrame()->GetWidgetName());
	//app->Script("pack %s -side bottom -anchor nw -expand n -padx 10 -pady 10", 
	//	m_pbStartMABMIS->GetWidgetName());


	/////////////////////////////////////////////////////////////////////


	// set file browser widget
	if (!m_fbwChooseDirectory)
		m_fbwChooseDirectory = vtkKWFileBrowserWidget::New();
	m_fbwChooseDirectory->SetParent(m_flData->GetFrame());
	m_fbwChooseDirectory->Create();
	m_fbwChooseDirectory->SetBorderWidth(2);
	m_fbwChooseDirectory->SetReliefToGroove();
	m_fbwChooseDirectory->SetPadX(2);
	m_fbwChooseDirectory->SetPadY(2);
	m_fbwChooseDirectory->SetWidth(540);
	m_fbwChooseDirectory->SetHeight(170);
	m_fbwChooseDirectory->SetFavoriteDirectoriesFrameVisibility(0);
	m_fbwChooseDirectory->SetFileListTableVisibility(1);
	m_fbwChooseDirectory->SetFocusToDirectoryExplorer();
	m_fbwChooseDirectory->MultipleSelectionOff();
	m_fbwChooseDirectory->GetFileListTable();
	//m_fbwChooseDirectory->GetDirectoryExplorer()->SetWidth(150);
	//m_fbwChooseDirectory->GetFileListTable()->SetWidth(250);
	//m_fbwChooseDirectory->FilterFilesByExtensions(".hdr .img");
	//m_fbwChooseDirectory->FilterFilesByExtensions(".txt .hdr .img");
	m_fbwChooseDirectory->FilterFilesByExtensions(".txt");
	m_fbwChooseDirectory->OpenDirectory(m_cCurrentDirectory);
	app->Script("pack %s -side top -anchor nw -expand 1 -fill both -padx 2 -pady 2 -in %s", 
		m_fbwChooseDirectory->GetWidgetName(), m_flData->GetFrame()->GetWidgetName());

	return 1;
}

void vtkKWMABMISGUI::ShowHelpAboutDialog()
{
	vtkKWMessageDialog * msg_dlg = vtkKWMessageDialog::New();
	msg_dlg->SetParent(this->window);
	msg_dlg->SetApplication(this->GetApplication());
	msg_dlg->Create();
	msg_dlg->SetText(" MABMIS GUI 1.0\n Author: Hongjun Jia (jiahj@med.unc.edu)\n Supported by: NIH grant R01EB006733 \n PI: Dinggang Shen \n Copyright (c) 2009-2010.");
	msg_dlg->Invoke();
	msg_dlg->Delete();
}

void vtkKWMABMISGUI::InitializePara()
{
	m_bFilesValid = false;
	m_cCurrentDirectory = new char[1024];
	return;
}
void vtkKWMABMISGUI::CheckFileValidity()
{
	m_bFilesValid = true;
	return;
}
void myitoa(int num, char* str, int digit)
{
	// translate the last (1~3) digits of a number to string
	if (digit == 1)
	{
		str[0] = num%10 + 48;
		str[1] = 0;
	}
	else if (digit == 2)
	{
		str[0] = (num%100)/10 + 48;
		str[1] = num%10 + 48;
		str[2] = 0;
	}
	else if (digit == 3)
	{
		str[0] = (num%1000)/100 + 48;
		str[1] = (num%100)/10 + 48;
		str[2] = num%10 + 48;
		str[3] = 0;
	}
	else
	{
		str[0] = 0;
	}
	return;
}

void vtkKWMABMISGUI::MABMISRunCallback()
{
	int atlas_size = 0;
	atlas_size = m_elAtlasesSize->GetWidget()->GetValueAsInt();

	// create necessary files name
	char MABMISFORGUINAME[MAX_FILE_NAME_LENGTH];
	MABMISFORGUINAME[MAX_FILE_NAME_LENGTH-1] = 0;
	strcpy(MABMISFORGUINAME, programFolder);

	char dataFolder[MAX_FILE_NAME_LENGTH];

#if defined(_WIN32) && !defined(__CYGWIN__)	// win32
	strcat(MABMISFORGUINAME, "/MABMIS_command.exe");
#else //unix64
	strcat(MABMISFORGUINAME, "/MABMIS_command");
#endif

	bool sanityCheckPass = true;
	int returnFromMABMIS = 0;

	int filelistsnum = m_fbwChooseDirectory->GetFileListTable()->GetNumberOfSelectedFileNames();

	if (filelistsnum == 0)
	{
		m_lWarningInfo->SetText("Please select one file with image IDs.");
		m_lWarningInfo->SetForegroundColor(1,0,0);
		m_lWarningInfo->SetFont("system 12 bold");
		return;
	} 
	else if (filelistsnum > 1)
	{
		m_lWarningInfo->SetText("Please select one file with image IDs.");
		m_lWarningInfo->SetForegroundColor(1,0,0);
		m_lWarningInfo->SetFont("system 12 bold");
		return;
	}
	else
	{		
		m_lWarningInfo->SetText("Select image ID lists file (.txt)!");
		m_lWarningInfo->SetForegroundColor(0,0,0);
		m_lWarningInfo->SetFont("system 12 normal");
		this->GetApplication()->ProcessPendingEvents();
	}

	char**	filesname = new char*[filelistsnum];
	for (int i = 0; i < filelistsnum; i++)
	{
		filesname[i] = new char[MAX_FILE_NAME_LENGTH];
		const char * tempfilename;
		tempfilename = m_fbwChooseDirectory->GetFileListTable()->GetNthSelectedFileName(i);
		strcpy(filesname[i], tempfilename);
	}
	
	strcpy(dataFolder, filesname[0]);
	for (int i = strlen(dataFolder); i>=0; i--)
	{
		if ((dataFolder[i]=='\\')||(dataFolder[i]=='/'))
		{
			dataFolder[i+1] = '\0';
			break;
		}
	}

	if (!vtksys::SystemTools::FileExists(MABMISFORGUINAME, true))
	{
#if defined(_WIN32) && !defined(__CYGWIN__)	// win32
		m_lWarningInfo->SetText("Missing executable file MABMIS_command.exe.");
#else //unix64
		m_lWarningInfo->SetText("Missing executable file MABMIS_command.");
#endif
		m_lWarningInfo->SetForegroundColor(1,0,0);
		m_lWarningInfo->SetFont("system 12 bold");
		//std::cerr << "Missing library file " << APLibName << "!" << std::endl;
		for (int i = 0; i < filelistsnum; i++) delete[] filesname[i];
		delete[] filesname;
		return;
	}

	// if there is space in the absolute path of data folder
	bool isSpaceInImageListFile = false;
	for (unsigned int i = 0; i < strlen(filesname[0]); i++)
	{
		if (filesname[0][i]==' ')
		{
			isSpaceInImageListFile = true;
			break;
		}
	}
	if (isSpaceInImageListFile)
	{
		m_lWarningInfo->SetText("Spaces found in the data file path!");
		m_lWarningInfo->SetForegroundColor(1,0,0);
		m_lWarningInfo->SetFont("system 12 bold");

		for (int i = 0; i < filelistsnum; i++) delete[] filesname[i];
		delete[] filesname;
		return;
	}

	if (!vtksys::SystemTools::FileExists(filesname[0], true))
	{
		m_lWarningInfo->SetText("Missing image ID file.");
		m_lWarningInfo->SetForegroundColor(1,0,0);
		m_lWarningInfo->SetFont("system 12 bold");
		//std::cerr << "No such file " << argv[1] << "!" << std::endl;
		for (int i = 0; i < filelistsnum; i++) delete[] filesname[i];
		delete[] filesname;
		return;
	}

	// load file numbers
	// 
	//std::cerr << "Loading file names and checking ... ";
	int filenumber = 0;
	std::ifstream fileNamesFile;
	fileNamesFile.open (filesname[0]);
	fileNamesFile.seekg (0, ios::end);// get length of file:
	filenumber = (unsigned int)(fileNamesFile.tellg()/3 + 0.5);
	fileNamesFile.seekg (0, ios::beg);
	fileNamesFile.close();
	fileNamesFile.open (filesname[0]);

	char** sub_ids = new char*[filenumber];
	for (int i=0; i< filenumber; i++)
	{
		//char temp_char;
		sub_ids[i] = new char[3];
		fileNamesFile >> sub_ids[i][0];
		fileNamesFile >> sub_ids[i][1];
		//idFile >> temp_char;
		sub_ids[i][2] = '\0';
	}
	fileNamesFile.close();

	// check if atlas_size is valid
	if ((atlas_size < 1)||(atlas_size > filenumber-1 ))
	{
		//
		m_lWarningInfo->SetText("The number of Atlases is invalid.");
		m_lWarningInfo->SetForegroundColor(1,0,0);
		m_lWarningInfo->SetFont("system 12 bold");
		//std::cerr << "No such file " << argv[1] << "!" << std::endl;
		for (int i = 0; i < filelistsnum; i++) delete[] filesname[i];
		delete[] filesname;
		for (int i = 0 ; i< filenumber; i++) delete[]  sub_ids[i];
		delete[] sub_ids;
		return;
	}

	//load file names and do sanity check

	char** imgfilenames = new char*[filenumber];
	char** segfilenames = new char*[atlas_size];
	for (int i = 0; i < filenumber; i++)
	{
		imgfilenames[i] = new char[MAX_FILE_NAME_LENGTH];
		
		strcpy(imgfilenames[i], dataFolder);
		strcat(imgfilenames[i], sub_ids[i]);
		strcat(imgfilenames[i], "_cbq_000.hdr");
	}
	for (int i = 0; i < atlas_size; i++)
	{
		segfilenames[i] = new char[MAX_FILE_NAME_LENGTH];

		strcpy(segfilenames[i], dataFolder);
		strcat(segfilenames[i], sub_ids[i]);
		strcat(segfilenames[i], "_seg_000.hdr");
	}


	// do sanity check 
	for (int i = 0; i < filenumber; i++)
	{
		char* tempfile1 = new char[1024];
		strcpy(tempfile1, imgfilenames[i]);
		int len1 = strlen(tempfile1);
		tempfile1[len1-1] = 'g';
		tempfile1[len1-2] = 'm';
		tempfile1[len1-3] = 'i';
		if (!vtksys::SystemTools::FileExists(imgfilenames[i], true))
		{
			sanityCheckPass = false;
			break;
		}
		if (!vtksys::SystemTools::FileExists(tempfile1, true))
		{
			sanityCheckPass = false;
			break;
		}
		delete[] tempfile1;
	}
	for (int i = 0; i < atlas_size; i++)
	{
		char* tempfile2 = new char[1024];
		strcpy(tempfile2, segfilenames[i]);
		int len2 = strlen(tempfile2);
		tempfile2[len2-1] = 'g';
		tempfile2[len2-2] = 'm';
		tempfile2[len2-3] = 'i';
		if (!vtksys::SystemTools::FileExists(segfilenames[i], true))
		{
			sanityCheckPass = false;
			break;
		}
		if (!vtksys::SystemTools::FileExists(tempfile2, true))
		{
			sanityCheckPass = false;
			break;
		}
		delete[] tempfile2;
	}

	if (!sanityCheckPass)
	{
		std::cerr << "Please check if all .hdr and .img file exist!" << std::endl;

		m_lWarningInfo->SetText("Please check all .hdr and .img files (ID_cbq_000)!");
		m_lWarningInfo->SetForegroundColor(1,0,0);
		m_lWarningInfo->SetFont("system 12 bold");

		for (int i = 0; i < filelistsnum; i++) delete[] filesname[i];
		delete[] filesname;
		for (int i = 0 ; i< filenumber; i++) delete[]  sub_ids[i];
		delete[] sub_ids;
		for (int i = 0 ; i< filenumber; i++) delete[]  imgfilenames[i];
		delete[] imgfilenames;
		for (int i = 0 ; i< atlas_size; i++) delete[]  segfilenames[i];
		delete[] segfilenames;

		return;
	}
	else
	{
		//std::cerr << "All paired .hdr .img files exist" << std::endl;
		m_lWarningInfo->SetText("MABMIS is running ...");
		m_lWarningInfo->SetForegroundColor(1,0,0);
		m_lWarningInfo->SetFont("system 12 bold");
		this->GetApplication()->ProcessPendingEvents();

		// start to process all valid inputs
	
		std::ostringstream oss (std::ostringstream::out);


		oss << "\"" <<MABMISFORGUINAME << "\""  << " "
			<< filesname[0] << " "
			<< atlas_size
			<< "	" << std::endl;

		std::cerr << oss.str().c_str();
		returnFromMABMIS = system(oss.str().c_str());

		if (returnFromMABMIS)  // EXIT_FAILURE == 1 
		{
			m_lWarningInfo->SetText("MABMIS failed!");
			m_lWarningInfo->SetForegroundColor(1,0,0);
			m_lWarningInfo->SetFont("system 12 bold");

			for (int i = 0; i < filelistsnum; i++) delete[] filesname[i];
			delete[] filesname;
			for (int i = 0 ; i< filenumber; i++) delete[]  sub_ids[i];
			delete[] sub_ids;
			for (int i = 0 ; i< filenumber; i++) delete[]  imgfilenames[i];
			delete[] imgfilenames;
			for (int i = 0 ; i< atlas_size; i++) delete[]  segfilenames[i];
			delete[] segfilenames;
			return;
		} 
		else  // EXIT_SUCCESS == 0
		{
			m_lWarningInfo->SetText("MABMIS is done successfully!");
			m_lWarningInfo->SetForegroundColor(0,0,0);
			m_lWarningInfo->SetFont("system 12 normal");
		}

	}

	// delete
	for (int i = 0; i < filelistsnum; i++) delete[] filesname[i];
	delete[] filesname;
	for (int i = 0 ; i< filenumber; i++) delete[]  sub_ids[i];
	delete[] sub_ids;
	for (int i = 0 ; i< filenumber; i++) delete[]  imgfilenames[i];
	delete[] imgfilenames;
	for (int i = 0 ; i< atlas_size; i++) delete[]  segfilenames[i];
	delete[] segfilenames;

	return;
}