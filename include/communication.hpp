//
// Created by xetql on 20.12.18.
//

#ifndef CA_ROAD_COMMUNICATION_HPP
#define CA_ROAD_COMMUNICATION_HPP

#include <mpi.h>
struct CommunicationDatatype {

    MPI_Datatype vec_datatype;
    MPI_Datatype elements_datatype;

    CommunicationDatatype(const MPI_Datatype &vec,
                          const MPI_Datatype &elements) :
        vec_datatype(vec), elements_datatype(elements){}
    void free_datatypes(){
        MPI_Type_free(&vec_datatype);
        MPI_Type_free(&elements_datatype);
    }
};
#endif //CA_ROAD_COMMUNICATION_HPP
