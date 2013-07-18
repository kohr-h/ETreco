/*  Copyright (c) 2007 Mark VanderVoord
 * 
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 * 
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.+
 * 
 * The end-user documentation included with the redistribution, if 
 * any, must include the following acknowledgment: "This product 
 * includes software developed for the CEXCeption Project, by Mark 
 * VanderVoord and other contributors", in the same place and form 
 * as other third-party acknowledgments.  Alternately, this 
 * acknowledgment may appear in the software itself, in the same 
 * form and location as other such third-party acknowledgments.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

#include <stdio.h>
#include "CException.h"

volatile CEXCEPTION_FRAME_T CExceptionFrames[CEXCEPTION_NUM_ID] = { {0} };

//------------------------------------------------------------------------------------------
//  Throw
//------------------------------------------------------------------------------------------

void Throw(CEXCEPTION_T ExceptionID)
{
    unsigned int MY_ID = CEXCEPTION_GET_ID;
    CExceptionFrames[MY_ID].Exception = ExceptionID;
    if (CExceptionFrames[MY_ID].pFrame)
    {
        longjmp(*CExceptionFrames[MY_ID].pFrame, 1);
    }
    CEXCEPTION_NO_CATCH_HANDLER(MY_ID);
}

