#ifndef UCV_CREATENEWKEY_h
#define UCV_CREATENEWKEY_h


#include <vtkm/worklet/WorkletMapField.h>
#include <vtkm/worklet/Keys.h>

struct CreateNewKeyWorklet : public vtkm::worklet::WorkletMapField
{

    VTKM_CONT
    CreateNewKeyWorklet(vtkm::Id rawDimx, vtkm::Id rawDimy, vtkm::Id rawDimz,
                        vtkm::Id numberBlockx, vtkm::Id numberBlocky, vtkm::Id numberBlockz,
                        vtkm::Id blocksize) : m_rawDimx(rawDimx), m_rawDimy(rawDimy), m_rawDimz(rawDimz),
                                              m_numberBlockx(numberBlockx), m_numberBlocky(numberBlocky), m_numberBlockz(numberBlockz),
                                              m_blocksize(blocksize){};

    vtkm::Id m_rawDimx;
    vtkm::Id m_rawDimy;
    vtkm::Id m_rawDimz;

    vtkm::Id m_numberBlockx;
    vtkm::Id m_numberBlocky;
    vtkm::Id m_numberBlockz;

    vtkm::Id m_blocksize;

    typedef void ControlSignature(FieldIn, FieldOut);
    typedef void ExecutionSignature(_1, _2);

    template <typename T>
    VTKM_EXEC void operator()(const T &globalIndex, vtkm::Id &hixelId) const
    {
        vtkm::Id globalId = static_cast<vtkm::Id>(globalIndex);
        // change globalId to the hixel index
        // compute the raw index firstly
        vtkm::Id rawx = globalId % m_rawDimx;
        vtkm::Id rawy = (globalId / m_rawDimx) % m_rawDimy;
        vtkm::Id rawz = globalId / m_rawDimx / m_rawDimy;

        vtkm::Id newidx = rawx / m_blocksize;
        vtkm::Id newidy = rawy / m_blocksize;
        vtkm::Id newidz = rawz / m_blocksize;

        hixelId = newidx + newidy * m_numberBlockx + newidz * m_numberBlockx * m_numberBlocky;
    };
};


#endif //UCV_CREATENEWKEY_h