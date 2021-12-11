#if defined(_MSC_VER)
#define NOMINMAX
#pragma once
#endif

#ifndef LOADOBJ_H
#define LOADOBJ_H

#include <string>
#include <fstream>
#include <sstream>
#include <iostream>

#include <vector>

#include <assimp/Importer.hpp>
#include <assimp/scene.h>
#include <assimp/postprocess.h>

#include "vec3.h"


//unsigned int TextureFromFile(const char* path, const std::string& directory, bool gamma = false);

class Mesh {
public:
    Mesh(int nv,int nf):v_num(nv),f_num(nf){
        v_pos = new Point3f[nv];
        f_indics = new int[3 * nf];
    }

    int v_num = 0;
    int f_num = 0;
    Point3f* v_pos = nullptr;
    int* f_indics = nullptr;
    Normal3f* vn = nullptr;
    Vector3f* vt = nullptr;
    Point2f* uv = nullptr;
};

class Model
{
public:
    // model data 
    //std::vector<Texture> textures_loaded;	// stores all the textures loaded so far, optimization to make sure textures aren't loaded more than once.
    std::vector<Mesh>    meshes;
    std::string directory;
    bool gammaCorrection;

    // constructor, expects a filepath to a 3D model.
    Model(std::string const& path, bool gamma = false) : gammaCorrection(gamma)
    {
        loadModel(path);
    }


private:

    // loads a model with supported ASSIMP extensions from file and stores the resulting meshes in the meshes vector.
    void loadModel(std::string const& path)
    {
        // read file via ASSIMP
        Assimp::Importer importer;
        const aiScene* scene = importer.ReadFile(path, aiProcess_Triangulate | aiProcess_GenSmoothNormals  /*| aiProcess_FlipUVs*/ | aiProcess_CalcTangentSpace);
        // check for errors
        if (!scene || scene->mFlags & AI_SCENE_FLAGS_INCOMPLETE || !scene->mRootNode) // if is Not Zero
        {
            std::cerr << "ERROR::ASSIMP:: " << importer.GetErrorString() << std::endl;
            return;
        }
        // retrieve the directory path of the filepath
        directory = path.substr(0, path.find_last_of('/'));

        // process ASSIMP's root node recursively
        processNode(scene->mRootNode, scene);
    }

    // processes a node in a recursive fashion. Processes each individual mesh located at the node and repeats this process on its children nodes (if any).
    void processNode(aiNode* node, const aiScene* scene)
    {
        // process each mesh located at the current node
        for (unsigned int i = 0; i < node->mNumMeshes; i++)
        {
            // the node object only contains indices to index the actual objects in the scene. 
            // the scene contains all the data, node is just to keep stuff organized (like relations between nodes).
            aiMesh* mesh = scene->mMeshes[node->mMeshes[i]];
            meshes.push_back(processMesh(mesh, scene));
        }
        // after we've processed all of the meshes (if any) we then recursively process each of the children nodes
        for (unsigned int i = 0; i < node->mNumChildren; i++)
        {
            processNode(node->mChildren[i], scene);
        }

    }

    Mesh processMesh(aiMesh* mesh, const aiScene* scene)
    {
        // data to fill
        //std::vector<Texture> textures;
        Mesh out_mesh(mesh->mNumVertices, mesh->mNumFaces);
        if (mesh->HasNormals())
            out_mesh.vn=new Normal3f[mesh->mNumVertices];
        if (mesh->mTextureCoords[0])
        {
            out_mesh.uv = new Point2f[mesh->mNumVertices];
            out_mesh.vt = new Vector3f[mesh->mNumVertices];
        }

        // walk through each of the mesh's vertices
        for (unsigned int i = 0; i < mesh->mNumVertices; i++)
        {
            // positions
            out_mesh.v_pos[i] = Point3f(mesh->mVertices[i].x, mesh->mVertices[i].y, mesh->mVertices[i].z);

            // normals
            if (mesh->HasNormals())
            {
                out_mesh.vn[i] = Normal3f(mesh->mNormals[i].x, mesh->mNormals[i].y, mesh->mNormals[i].z);
            }

            // texture coordinates
            if (mesh->mTextureCoords[0]) // does the mesh contain texture coordinates?
            {
                out_mesh.uv[i] = Point2f(mesh->mTextureCoords[0][i].x, mesh->mTextureCoords[0][i].y);
                // tangent
                out_mesh.vt[i] = Vector3f(mesh->mTangents[i].x, mesh->mTangents[i].y, mesh->mTangents[i].z);
                // a vertex can contain up to 8 different texture coordinates. We thus make the assumption that we won't 
                // use models where a vertex can have multiple texture coordinates so we always take the first set (0).
            }
            //else
                //vertex.TexCoords = glm::vec2(0.0f, 0.0f);
        }
        // now wak through each of the mesh's faces (a face is a mesh its triangle) and retrieve the corresponding vertex indices.
        for (unsigned int i = 0; i < mesh->mNumFaces; i++)
        {
            aiFace face = mesh->mFaces[i];
            // retrieve all indices of the face and store them in the indices vector
            for (unsigned int j = 0; j < face.mNumIndices; j++)
                out_mesh.f_indics[3 * i + j] = face.mIndices[j];
        }
        // process materials
        aiMaterial* material = scene->mMaterials[mesh->mMaterialIndex];
        // we assume a convention for sampler names in the shaders. Each diffuse texture should be named
        // as 'texture_diffuseN' where N is a sequential number ranging from 1 to MAX_SAMPLER_NUMBER. 
        // Same applies to other texture as the following list summarizes:
        // diffuse: texture_diffuseN
        // specular: texture_specularN
        // normal: texture_normalN

        //// 1. diffuse maps
        //vector<Texture> diffuseMaps = loadMaterialTextures(material, aiTextureType_DIFFUSE, "texture_diffuse");
        //textures.insert(textures.end(), diffuseMaps.begin(), diffuseMaps.end());
        //// 2. specular maps
        //vector<Texture> specularMaps = loadMaterialTextures(material, aiTextureType_SPECULAR, "texture_specular");
        //textures.insert(textures.end(), specularMaps.begin(), specularMaps.end());
        //// 3. normal maps
        //std::vector<Texture> normalMaps = loadMaterialTextures(material, aiTextureType_HEIGHT, "texture_normal");
        //textures.insert(textures.end(), normalMaps.begin(), normalMaps.end());
        //// 4. height maps
        //std::vector<Texture> heightMaps = loadMaterialTextures(material, aiTextureType_AMBIENT, "texture_height");
        //textures.insert(textures.end(), heightMaps.begin(), heightMaps.end());

        // return a mesh object created from the extracted mesh data
        return out_mesh;
        //return Mesh(vertices, indices, textures);

    }

    // checks all material textures of a given type and loads the textures if they're not loaded yet.
    // the required info is returned as a Texture struct.
    //vector<Texture> loadMaterialTextures(aiMaterial* mat, aiTextureType type, string typeName)
    //{
    //    vector<Texture> textures;
    //    for (unsigned int i = 0; i < mat->GetTextureCount(type); i++)
    //    {
    //        aiString str;
    //        mat->GetTexture(type, i, &str);
    //        // check if texture was loaded before and if so, continue to next iteration: skip loading a new texture
    //        bool skip = false;
    //        for (unsigned int j = 0; j < textures_loaded.size(); j++)
    //        {
    //            if (std::strcmp(textures_loaded[j].path.data(), str.C_Str()) == 0)
    //            {
    //                textures.push_back(textures_loaded[j]);
    //                skip = true; // a texture with the same filepath has already been loaded, continue to next one. (optimization)
    //                break;
    //            }
    //        }
    //        if (!skip)
    //        {   // if texture hasn't been loaded already, load it
    //            Texture texture;
    //            texture.id = TextureFromFile(str.C_Str(), this->directory);
    //            texture.type = typeName;
    //            texture.path = str.C_Str();
    //            textures.push_back(texture);
    //            textures_loaded.push_back(texture);  // store it as texture loaded for entire model, to ensure we won't unnecesery load duplicate textures.
    //        }
    //    }
    //    return textures;
    //}
};


//unsigned int TextureFromFile(const char* path, const string& directory, bool gamma)
//{
//    string filename = string(path);
//    filename = directory + '/' + filename;
//
//    unsigned int textureID;
//    glGenTextures(1, &textureID);
//
//    int width, height, nrComponents;
//    unsigned char* data = stbi_load(filename.c_str(), &width, &height, &nrComponents, 0);
//    if (data)
//    {
//        GLenum format;
//        if (nrComponents == 1)
//            format = GL_RED;
//        else if (nrComponents == 3)
//            format = GL_RGB;
//        else if (nrComponents == 4)
//            format = GL_RGBA;
//
//        glBindTexture(GL_TEXTURE_2D, textureID);
//        glTexImage2D(GL_TEXTURE_2D, 0, format, width, height, 0, format, GL_UNSIGNED_BYTE, data);
//        glGenerateMipmap(GL_TEXTURE_2D);
//
//        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
//        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
//        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
//        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
//
//        stbi_image_free(data);
//    }
//    else
//    {
//        std::cout << "Texture failed to load at path: " << path << std::endl;
//        stbi_image_free(data);
//    }
//
//    return textureID;
//}

#endif // !LOADOBJ_H

