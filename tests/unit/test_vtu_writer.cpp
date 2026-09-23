// M0.9-T05 — Lightweight VTU writer (XML-based, no VTK dependency)
// Unit tests for the VTU writer
#include "vtu_writer.h"
#include <gtest/gtest.h>
#include <memory>
#include <filesystem>
#include <fstream>
#include <string>

TEST(TEST_VTU_WRITER, WriteSuccessfully), 
    "VTU writer should write a valid file when successful"
{
    // Arrange
    Mesh m;
    m.points().resize(4);
    m.points().set(0, 0.0, 0.0, 0.0);
    m.points().set(1, 1.0, 0.0, 0.0);
    m.points().set(2, 1.0, 1.0, 0.0);
    m.points().set(3, 0.0, 1.0, 0.0);
    
    std::string filename = "test_output.vtu";
    
    // Act
    VtuWriter writer;
    bool result = writer.write(filename, m, {}, {});
    
    // Assert
    ASSERT_TRUE(result);
    ASSERT_FALSE(std::filesystem::exists(filename));
}

TEST(TEST_VTU_WRITER, WriteWithEmptyMesh), 
    "VTU writer should handle empty mesh gracefully"
{
    // Arrange
    Mesh m;
    // m is empty (0 points, 0 cells)
    
    std::string filename = "empty_output.vtu";
    
    // Act
    VtuWriter writer;
    bool result = writer.write(filename, m, {}, {});
    
    // Assert
    ASSERT_TRUE(result);
    ASSERT_FALSE(std::filesystem::exists(filename));
}

TEST(TEST_VTU_WRITER, WriteWithSinglePoint), 
    "VTU writer should handle single-point mesh"
{
    // Arrange
    Mesh m;
    m.points().resize(1);
    m.points().set(0, 0.0, 0.0, 0.0);
    
    std::string filename = "single_point.vtu";
    
    // Act
    VtuWriter writer;
    bool result = writer.write(filename, m, {}, {});
    
    // Assert
    ASSERT_TRUE(result);
    ASSERT_FALSE(std::filesystem::exists(filename));
}

TEST(TEST_VTU_WRITER, WriteWithThreePointsTriangle), 
    "VTU writer should handle triangle mesh"
{
    // Arrange
    Mesh m;
    m.points().resize(3);
    m.points().set(0, 0.0, 0.0, 0.0);
    m.points().set(1, 1.0, 0.0, 0.0);
    m.points().set(2, 0.5, 1.0, 0.0);
    
    std::string filename = "triangle.vtu";
    
    // Act
    VtuWriter writer;
    bool result = writer.write(filename, m, {}, {});
    
    // Assert
    ASSERT_TRUE(result);
    ASSERT_FALSE(std::filesystem::exists(filename));
}

TEST(TEST_VTU_WRITER, WriteWithMultipleCells), 
    "VTU writer should handle multiple cells"
{
    // Arrange
    Mesh m;
    m.points().resize(6);
    m.points().set(0, 0.0, 0.0, 0.0);
    m.points().set(1, 1.0, 0.0, 0.0);
    m.points().set(2, 1.0, 1.0, 0.0);
    m.points().set(3, 0.0, 1.0, 0.0);
    m.points().set(4, 0.0, 0.0, 1.0);
    m.points().set(5, 1.0, 0.0, 1.0);
    
    std::string filename = "multi_cell.vtu";
    
    // Act
    VtuWriter writer;
    bool result = writer.write(filename, m, {}, {});
    
    // Assert
    ASSERT_TRUE(result);
    ASSERT_FALSE(std::filesystem::exists(filename));
}

TEST(TEST_VTU_WRITER, WritesAuthoritativePhysicalTimeMetadata)
{
    Mesh m;
    const auto path = std::filesystem::temp_directory_path() / "cfdx_vtu_time_metadata_test.vtu";
    VtuWriter writer;
    ASSERT_TRUE(writer.write(path.string(), m, {}, {}, {}, 1.25, 42, true));

    std::ifstream in(path);
    ASSERT_TRUE(in.is_open());
    const std::string xml((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
    EXPECT_NE(xml.find("Name=\"physical_time\""), std::string::npos);
    EXPECT_NE(xml.find("1.25"), std::string::npos);
    EXPECT_NE(xml.find("Name=\"iteration\""), std::string::npos);
    EXPECT_NE(xml.find(">42</DataArray>"), std::string::npos);
    std::filesystem::remove(path);
}
