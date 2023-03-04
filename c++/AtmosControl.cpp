#include "AtmosControl.h"

AtmosControl::AtmosControl()
{
    time = 0.0;
    precipFact = 0.5;
    petFact = 1.0;
}

double AtmosControl::get_precip()
{
    double timeStepData = 300.0;
    size_t timeStep = time / timeStepData;
    double dataVal = 0.0;
    size_t dataIndex = 0;
    if (timeStep >= 0 && timeStep < atmosForcingData.size())
    {
        // mm to m conversion done below.
        if (atmosForcingData.at(timeStep).at(dataIndex) > 0.0)
        {
            dataVal = 0.001 * precipFact * atmosForcingData.at(timeStep).at(dataIndex) / timeStepData;
        }
    }
    return dataVal;
}

double AtmosControl::get_pet()
{
    double timeStepData = 300.0;
    size_t timeStep = time / timeStepData;
    double dataVal = 0.0;
    size_t dataIndex = 1;
    if (timeStep >= 0 && timeStep < atmosForcingData.size())
    {
        // mm to m conversion done below.
        if (atmosForcingData.at(timeStep).at(dataIndex) > 0.0)
        {
            dataVal = 0.001 * petFact * atmosForcingData.at(timeStep).at(dataIndex) / timeStepData;
        }
    }
    return dataVal;
}

double AtmosControl::get_air_temp()
{
    double timeStepData = 300.0;
    size_t timeStep = time / timeStepData;
    double dataVal = 0.0;
    size_t dataIndex = 2;
    if (timeStep >= 0 && timeStep < atmosForcingData.size())
    {
        if (atmosForcingData.at(timeStep).at(dataIndex) > 0.0)
        {
            dataVal = atmosForcingData.at(timeStep).at(dataIndex);
        }
    }
    return dataVal;
}

void AtmosControl::set_atmos_forcing_data(std::vector < std::vector<std::string> > atmosForcingDataStr)
{
    // Save and remove the data header row.
    dataHeaders = atmosForcingDataStr.at(0);
    atmosForcingDataStr.erase(atmosForcingDataStr.begin());
    dataHeaders.erase(dataHeaders.begin()); // remove date time header
    // Convert the data.
    for (size_t i = 0; i < atmosForcingDataStr.size(); i++)
    {
        // Extract datetime column into a separate vector.
        dateTime.push_back(atmosForcingDataStr.at(i).at(0));
        atmosForcingDataStr.at(i).erase(atmosForcingDataStr.at(i).begin());
        // Convert data from string to double type.
        std::vector<double> dataRow;
        for (size_t j = 0; j < atmosForcingDataStr.at(i).size(); j++)
        {
            dataRow.push_back(atof(atmosForcingDataStr.at(i).at(j).c_str()));
        }
        atmosForcingData.push_back(dataRow);
    }
}

void AtmosControl::print_atmos_forcing_data()
{
    // Print the header.
    std::cout << "\nAtmospheric forcing data:";
    std::cout << "\n\ndatetime [datetime]";
    for (size_t i = 0; i < dataHeaders.size(); i++)
    {
        std::cout << dataHeaders.at(i);
        if (i < dataHeaders.size() - 1)
        {
            std::cout << "\t";
        }
    }
    // Print the data.
    for (size_t i = 0; i < atmosForcingData.size(); i++)
    {
        if ((i >= 0 && i < 5) || i > atmosForcingData.size() - 6)
        {
            std::cout << "\n";
            std::cout << dateTime.at(i);
            std::cout << "\t";
        }
        if (i < 5)
        {
            for (size_t j = 0; j < atmosForcingData.at(i).size(); j++)
            {
                std::cout << atmosForcingData.at(i).at(j);
                if (j < atmosForcingData.at(i).size() - 1)
                {
                    std::cout << "\t";
                }
            }
        }
        else if (i == 5)
        {
            std::cout << "\n";
            std::cout << "...";
        }
        else if (i > atmosForcingData.size() - 6)
        {
            for (size_t j = 0; j < atmosForcingData.at(i).size(); j++)
            {
                std::cout << atmosForcingData.at(i).at(j);
                if (j < atmosForcingData.at(i).size() - 1)
                {
                    std::cout << "\t";
                }
            }
        }
    }
    std::cout << "\n";
}